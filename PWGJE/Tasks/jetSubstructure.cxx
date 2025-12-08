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

// jet analysis tasks (subscribing to jet finder task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "PWGJE/DataModel/JetSubstructure.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include <math.h>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetSubstructureTask {
  Produces<aod::CJetSSs> jetSubstructureDataTable;
  Produces<aod::CMCDJetSSs> jetSubstructureMCDTable;
  Produces<aod::CMCPJetSSs> jetSubstructureMCPTable;
  Produces<aod::CEWSJetSSs> jetSubstructureDataSubTable;

  Produces<aod::ChargedSPs> jetSplittingsDataTable;
  Produces<aod::ChargedMCDetectorLevelSPs> jetSplittingsMCDTable;
  Produces<aod::ChargedMCParticleLevelSPs> jetSplittingsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedSPs> jetSplittingsDataSubTable;

  Produces<aod::ChargedPRs> jetPairsDataTable;
  Produces<aod::ChargedMCDetectorLevelPRs> jetPairsMCDTable;
  Produces<aod::ChargedMCParticleLevelPRs> jetPairsMCPTable;
  Produces<aod::ChargedEventWiseSubtractedPRs> jetPairsDataSubTable;

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};
  Configurable<bool> doPairBkg{"doPairBkg", true, "save bkg pairs"};
  Configurable<float> pairConstituentPtMin{"pairConstituentPtMin", 1.0, "pt cut off for constituents going into pairs"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in jet finding"};
  Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
  Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied to jet finding if applyTrackingEfficiency is true"};

  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  std::vector<float> energyMotherVec;
  std::vector<float> ptLeadingVec;
  std::vector<float> ptSubLeadingVec;
  std::vector<float> thetaVec;
  std::vector<float> nSub;
  std::vector<float> pairJetPtVec;
  std::vector<float> pairJetEnergyVec;
  std::vector<float> pairJetThetaVec;
  std::vector<float> pairJetPerpCone1PtVec;
  std::vector<float> pairJetPerpCone1EnergyVec;
  std::vector<float> pairJetPerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone1PtVec;
  std::vector<float> pairPerpCone1PerpCone1EnergyVec;
  std::vector<float> pairPerpCone1PerpCone1ThetaVec;
  std::vector<float> pairPerpCone1PerpCone2PtVec;
  std::vector<float> pairPerpCone1PerpCone2EnergyVec;
  std::vector<float> pairPerpCone1PerpCone2ThetaVec;
  float angularity;
  float leadingConstituentPt;
  float perpConeRho;

  HistogramRegistry registry;

  int trackSelection = -1;

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
    jetReclusterer.ghostRepeatN = 0;

    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (applyTrackingEfficiency) {
      if (trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (trackingEfficiency->size() + 1 != trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "jetFinder workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  Preslice<aod::JetTracks> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::JetTracksSub> TracksPerCollisionDataSub = aod::bkgcharged::collisionId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  template <bool isMCP, bool isSubtracted, typename T, typename U>
  void jetReclustering(T const& jet, U& splittingTable)
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

    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      std::vector<int32_t> tracks;
      std::vector<int32_t> candidates;
      std::vector<int32_t> clusters;
      for (const auto& constituent : sorted_by_pt(parentSubJet2.constituents())) {
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          tracks.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
      }
      splittingTable(jet.globalIndex(), tracks, clusters, candidates, parentSubJet2.perp(), parentSubJet2.eta(), parentSubJet2.phi(), 0);
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      energyMotherVec.push_back(daughterSubJet.e());
      ptLeadingVec.push_back(parentSubJet1.pt());
      ptSubLeadingVec.push_back(parentSubJet2.pt());
      thetaVec.push_back(theta);

      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          auto zg = z;
          auto rg = theta;
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

  template <bool isMC, typename T, typename U, typename V, typename M>
  void jetPairing(T const& jet, U const& tracks, V const& slicer, M& pairTable)
  {
    pairJetPtVec.clear();
    pairJetEnergyVec.clear();
    pairJetThetaVec.clear();
    std::vector<typename U::iterator> tracksVec;
    std::vector<int32_t> tracksVecIds;
    for (auto const& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= pairConstituentPtMin) {
        tracksVec.push_back(constituent);
        tracksVecIds.push_back(constituent.globalIndex());
      }
    }
    if (tracksVec.size() >= 1) {
      for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
        for (typename std::vector<typename U::iterator>::size_type track2Index = track1Index + 1; track2Index < tracksVec.size(); track2Index++) {
          pairJetPtVec.push_back(tracksVec.at(track1Index).pt() * tracksVec.at(track2Index).pt());
          pairJetEnergyVec.push_back(2.0 * tracksVec.at(track1Index).energy() * tracksVec.at(track2Index).energy());
          pairJetThetaVec.push_back(jetutilities::deltaR(tracksVec.at(track1Index), tracksVec.at(track2Index)));
          pairTable(jet.globalIndex(), tracksVecIds.at(track1Index), tracksVecIds.at(track2Index), -1, -1);
        }
      }
    }

    pairJetPerpCone1PtVec.clear();
    pairJetPerpCone1EnergyVec.clear();
    pairJetPerpCone1ThetaVec.clear();
    pairPerpCone1PerpCone1PtVec.clear();
    pairPerpCone1PerpCone1EnergyVec.clear();
    pairPerpCone1PerpCone1ThetaVec.clear();
    pairPerpCone1PerpCone2PtVec.clear();
    pairPerpCone1PerpCone2EnergyVec.clear();
    pairPerpCone1PerpCone2ThetaVec.clear();

    int32_t collisionId = -1;
    if constexpr (!isMC) {
      collisionId = jet.collisionId();
    } else {
      collisionId = jet.mcCollisionId();
    }
    auto tracksPerCollision = tracks.sliceBy(slicer, collisionId);

    float perpCone1Phi = RecoDecay::constrainAngle<float, float>(jet.phi() + (M_PI / 2.));
    float perpCone2Phi = RecoDecay::constrainAngle<float, float>(jet.phi() - (M_PI / 2.));
    float perpCone1Pt = 0.0;
    float perpCone2Pt = 0.0;
    std::vector<typename U::iterator> tracksPerpCone1Vec;
    std::vector<typename U::iterator> tracksPerpCone2Vec;
    for (auto const& track : tracksPerCollision) {
      if (!jetderiveddatautilities::applyTrackKinematics(track)) {
        continue;
      }

      if constexpr (!std::is_same_v<std::decay_t<U>, aod::JetParticles>) {
        if (!jetfindingutilities::isTrackSelected<typename U::iterator, typename U::iterator>(track, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning)) {
          continue;
        }
      }

      float deltaPhi1 = track.phi() - perpCone1Phi;
      deltaPhi1 = RecoDecay::constrainAngle<float, float>(deltaPhi1, -M_PI);
      float deltaPhi2 = track.phi() - perpCone2Phi;
      deltaPhi2 = RecoDecay::constrainAngle<float, float>(deltaPhi2, -M_PI);
      float deltaEta = jet.eta() - track.eta();

      if (TMath::Sqrt((deltaPhi1 * deltaPhi1) + (deltaEta * deltaEta)) <= jet.r() / 100.0) {
        if (track.pt() >= pairConstituentPtMin) {
          tracksPerpCone1Vec.push_back(track);
        }
        perpCone1Pt += track.pt();
      }
      if (TMath::Sqrt((deltaPhi2 * deltaPhi2) + (deltaEta * deltaEta)) <= jet.r() / 100.0) {
        if (track.pt() >= pairConstituentPtMin) {
          tracksPerpCone2Vec.push_back(track);
        }
        perpCone2Pt += track.pt();
      }
    }
    perpConeRho = (perpCone1Pt + perpCone2Pt) / (2 * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)); // currently done per jet - could be better to do for leading jet if pushing to very low pT
    if (doPairBkg) {
      if (tracksVec.size() >= 1 && tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairJetPerpCone1PtVec.push_back(tracksVec.at(track1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            pairJetPerpCone1EnergyVec.push_back(2.0 * tracksVec.at(track1Index).energy() * tracksPerpCone1Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle(tracksVec.at(track1Index).phi() - (tracksPerpCone1Vec.at(track2Index).phi() - (M_PI / 2.)), -M_PI);
            float dEta = tracksVec.at(track1Index).eta() - tracksPerpCone1Vec.at(track2Index).eta();
            pairJetPerpCone1ThetaVec.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
        }
      }

      if (tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksPerpCone1Vec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = track1Index + 1; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairPerpCone1PerpCone1PtVec.push_back(tracksPerpCone1Vec.at(track1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            pairPerpCone1PerpCone1EnergyVec.push_back(2.0 * tracksPerpCone1Vec.at(track1Index).energy() * tracksPerpCone1Vec.at(track2Index).energy());
            pairPerpCone1PerpCone1ThetaVec.push_back(jetutilities::deltaR(tracksPerpCone1Vec.at(track1Index), tracksPerpCone1Vec.at(track2Index)));
          }
        }
      }

      if (tracksPerpCone1Vec.size() >= 1 && tracksPerpCone2Vec.size() >= 1) {
        for (typename std::vector<typename U::iterator>::size_type track1Index = 0; track1Index < tracksPerpCone1Vec.size(); track1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone2Vec.size(); track2Index++) {
            pairPerpCone1PerpCone2PtVec.push_back(tracksPerpCone1Vec.at(track1Index).pt() * tracksPerpCone2Vec.at(track2Index).pt());
            pairPerpCone1PerpCone2EnergyVec.push_back(2.0 * tracksPerpCone1Vec.at(track1Index).energy() * tracksPerpCone2Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle((tracksPerpCone1Vec.at(track1Index).phi() - (M_PI / 2.)) - (tracksPerpCone2Vec.at(track2Index).phi() + (M_PI / 2.)), -M_PI);
            float dEta = tracksPerpCone1Vec.at(track1Index).eta() - tracksPerpCone2Vec.at(track2Index).eta();
            pairPerpCone1PerpCone2ThetaVec.push_back(std::sqrt(dEta * dEta + dPhi * dPhi));
          }
        }
      }
    }
  }

  template <typename T, typename U>
  void jetSubstructureSimple(T const& jet, U const& /*tracks*/)
  {
    angularity = 0.0;
    leadingConstituentPt = 0.0;
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= leadingConstituentPt) {
        leadingConstituentPt = constituent.pt();
      }
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent), alpha);
    }
    angularity /= (jet.pt() * (jet.r() / 100.f));
  }

  template <bool isSubtracted, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& jet, U const& tracks, V const& trackSlicer, M& outputTable, N& splittingTable, O& pairTable)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, tracks, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet, splittingTable);
    jetPairing<false>(jet, tracks, trackSlicer, pairTable);
    jetSubstructureSimple(jet, tracks);
    outputTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetSubstructureTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
                              aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable, jetPairsDataTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>::iterator const& jet,
                                          aod::JetTracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, TracksPerCollisionDataSub, jetSubstructureDataSubTable, jetSplittingsDataSubTable, jetPairsDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(typename soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, TracksPerCollision, jetSubstructureMCDTable, jetSplittingsMCDTable, jetPairsMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCD, "charged jet substructure", false);

  void processChargedJetsMCP(typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                             aod::JetParticles const& particles)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, particles, particles, particles, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<true, false>(jet, jetSplittingsMCPTable);
    jetPairing<true>(jet, particles, ParticlesPerMcCollision, jetPairsMCPTable);
    jetSubstructureSimple(jet, particles);
    jetSubstructureMCPTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure"})};
}
