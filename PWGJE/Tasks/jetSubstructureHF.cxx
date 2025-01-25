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

// heavy-flavour jet substructure task (subscribing to jet finder hf task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTableData, typename JetTableMCD, typename JetTableMCP, typename JetTableDataSub, typename CandidateTable, typename CandidateTableMCP, typename SubstructureTableData, typename SubstructureTableMCD, typename SubstructureTableMCP, typename SubstructureTableDataSub, typename TracksSub>
struct JetSubstructureHFTask {
  Produces<SubstructureTableData> jetSubstructureDataTable;
  Produces<SubstructureTableMCD> jetSubstructureMCDTable;
  Produces<SubstructureTableMCP> jetSubstructureMCPTable;
  Produces<SubstructureTableDataSub> jetSubstructureDataSubTable;

  // Jet level configurables
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};
  Configurable<float> pairConstituentPtMin{"pairConstituentPtMin", 1.0, "pt cut off for constituents going into pairs"};

  Service<o2::framework::O2DatabasePDG> pdg;
  float candMass;

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<float> thetaVec;
  std::vector<float> nSub;
  std::vector<float> pairPtVec;
  std::vector<float> pairEnergyVec;
  std::vector<float> pairThetaVec;
  float angularity;
  float leadingConstituentPt;

  HistogramRegistry registry;
  void init(InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;

    candMass = jetcandidateutilities::getTablePDGMass<CandidateTable>();
  }

  template <bool isMCP, bool isSubtracted, typename T>
  void jetReclustering(T const& jet)
  {
    energyMotherVec.clear();
    ptLeadingVec.clear();
    ptSubLeadingVec.clear();
    thetaVec.clear();
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {

      bool isHFInSubjet1 = false;
      for (auto& subjet1Constituent : parentSubJet1.constituents()) {
        if (subjet1Constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidate)) {
          isHFInSubjet1 = true;
          break;
        }
      }
      if (!isHFInSubjet1) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      energyMotherVec.push_back(daughterSubJet.e());
      ptLeadingVec.push_back(parentSubJet1.pt());
      ptSubLeadingVec.push_back(parentSubJet2.pt());
      thetaVec.push_back(theta);
      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          if constexpr (!isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg);
          }
          if constexpr (!isSubtracted && isMCP) {
            registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg);
          }
          if constexpr (isSubtracted && !isMCP) {
            registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"), jet.pt(), zg);
            registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"), jet.pt(), rg);
          }
          softDropped = true;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
    }
    if constexpr (!isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd);
    }
  }

  template <typename T, typename U, typename V>
  void jetPairing(T const& jet, U const& /*tracks*/, V const& /*candidates*/)
  {
    pairPtVec.clear();
    pairEnergyVec.clear();
    pairThetaVec.clear();
    std::vector<std::decay_t<typename U::iterator>> tracksVec;
    std::vector<std::decay_t<typename V::iterator>> candidatesVec;
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= pairConstituentPtMin) {
        tracksVec.push_back(constituent);
      }
    }
    for (auto& candidate : jet.template candidates_as<V>()) {
      candidatesVec.push_back(candidate);
    }
    if (tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track2Index = 0; track2Index < tracksVec.size(); track2Index++) {
          pairPtVec.push_back(tracksVec.at(track1Index).pt() * tracksVec.at(track2Index).pt());
          pairEnergyVec.push_back(tracksVec.at(track1Index).energy() * tracksVec.at(track2Index).energy());
          pairThetaVec.push_back(jetutilities::deltaR(tracksVec.at(track1Index), tracksVec.at(track2Index)));
        }
      }
    }
    if (candidatesVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate1Index = 0; candidate1Index < candidatesVec.size(); candidate1Index++) {
        for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate2Index = 0; candidate2Index < candidatesVec.size(); candidate2Index++) {
          pairPtVec.push_back(candidatesVec.at(candidate1Index).pt() * candidatesVec.at(candidate2Index).pt());
          auto candidate1Energy = std::sqrt((candidatesVec.at(candidate1Index).p() * candidatesVec.at(candidate1Index).p()) + (candMass * candMass));
          auto candidate2Energy = std::sqrt((candidatesVec.at(candidate2Index).p() * candidatesVec.at(candidate2Index).p()) + (candMass * candMass));
          pairEnergyVec.push_back(candidate1Energy * candidate2Energy);
          pairThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidate1Index), candidatesVec.at(candidate2Index)));
        }
      }
    }
    if (candidatesVec.size() >= 1 && tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidateIndex = 0; candidateIndex < candidatesVec.size(); candidateIndex++) { // could just directly get the candidate and tracks here but keeping it consistent with above
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type trackIndex = 0; trackIndex < tracksVec.size(); trackIndex++) {
          pairPtVec.push_back(candidatesVec.at(candidateIndex).pt() * tracksVec.at(trackIndex).pt());
          auto candidateEnergy = std::sqrt((candidatesVec.at(candidateIndex).p() * candidatesVec.at(candidateIndex).p()) + (candMass * candMass));
          pairEnergyVec.push_back(candidateEnergy * tracksVec.at(trackIndex).energy());
          pairThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidateIndex), tracksVec.at(trackIndex)));
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void jetSubstructureSimple(T const& jet, U const& /*tracks*/, V const& /*candidates*/)
  {
    angularity = 0.0;
    leadingConstituentPt = 0.0;
    for (auto& candidate : jet.template candidates_as<V>()) {
      if (candidate.pt() >= leadingConstituentPt) {
        leadingConstituentPt = candidate.pt();
      }
      angularity += std::pow(candidate.pt(), kappa) * std::pow(jetutilities::deltaR(jet, candidate), alpha);
    }
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= leadingConstituentPt) {
        leadingConstituentPt = constituent.pt();
      }
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent), alpha);
    }
    angularity /= (jet.pt() * (jet.r() / 100.f));
  }

  template <bool isSubtracted, typename T, typename U, typename V, typename M>
  void analyseCharged(T const& jet, U const& tracks, V const& candidates, M& outputTable)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    for (auto& jetHFCandidate : jet.template candidates_as<V>()) { // should only be one at the moment
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidate), candMass);
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, candidates, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet);
    jetPairing(jet, tracks, candidates);
    jetSubstructureSimple(jet, tracks, candidates);
    outputTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairPtVec, pairEnergyVec, pairThetaVec, angularity, leadingConstituentPt);
  }

  void processChargedJetsData(typename JetTableData::iterator const& jet,
                              CandidateTable const& candidates,
                              aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, jetSubstructureDataTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsData, "HF jet substructure on data", false);

  void processChargedJetsDataSub(typename JetTableDataSub::iterator const& jet,
                                 CandidateTable const& candidates,
                                 TracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, candidates, jetSubstructureDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsDataSub, "HF jet substructure on data", false);

  void processChargedJetsMCD(typename JetTableMCD::iterator const& jet,
                             CandidateTable const& candidates,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, jetSubstructureMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCD, "HF jet substructure on data", false);

  void processChargedJetsMCP(typename JetTableMCP::iterator const& jet,
                             aod::JetParticles const& particles,
                             CandidateTableMCP const& candidates)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    for (auto& jetHFCandidate : jet.template candidates_as<CandidateTableMCP>()) {
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidate), candMass);
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, particles, particles, candidates, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<true, false>(jet);
    jetPairing(jet, particles, candidates);
    jetSubstructureSimple(jet, particles, candidates);
    jetSubstructureMCPTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairPtVec, pairEnergyVec, pairThetaVec, angularity, leadingConstituentPt);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCP, "HF jet substructure on MC particle level", false);
};
