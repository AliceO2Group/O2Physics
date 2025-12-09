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

#include "RecoDecay.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>

#include <TMath.h>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>

#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!

template <typename JetTableData, typename JetTableMCD, typename JetTableMCP, typename JetTableDataSub, typename CandidateTable, typename CandidateTableMCP, typename SubstructureTableData, typename SplittingsTableData, typename PairsTableData, typename SubstructureTableMCD, typename SplittingsTableMCD, typename PairsTableMCD, typename SubstructureTableMCP, typename SplittingsTableMCP, typename PairsTableMCP, typename SubstructureTableDataSub, typename SplittingsTableDataSub, typename PairsTableDataSub, typename TracksSub>
struct JetSubstructureHFTask {
  Produces<SubstructureTableData> jetSubstructureDataTable;
  Produces<SubstructureTableMCD> jetSubstructureMCDTable;
  Produces<SubstructureTableMCP> jetSubstructureMCPTable;
  Produces<SubstructureTableDataSub> jetSubstructureDataSubTable;

  Produces<SplittingsTableData> jetSplittingsDataTable;
  Produces<SplittingsTableMCD> jetSplittingsMCDTable;
  Produces<SplittingsTableMCP> jetSplittingsMCPTable;
  Produces<SplittingsTableDataSub> jetSplittingsDataSubTable;

  Produces<PairsTableData> jetPairsDataTable;
  Produces<PairsTableMCD> jetPairsMCDTable;
  Produces<PairsTableMCP> jetPairsMCPTable;
  Produces<PairsTableDataSub> jetPairsDataSubTable;

  // Jet level configurables
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
  float candMass;

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

    candMass = jetcandidateutilities::getTablePDGMass<CandidateTable>();

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
  PresliceOptional<aod::JTrackD0Subs> TracksPerD0DataSub = aod::bkgd0::candidateId;
  PresliceOptional<aod::JTrackDplusSubs> TracksPerDplusDataSub = aod::bkgdplus::candidateId;
  PresliceOptional<aod::JTrackDsSubs> TracksPerDsDataSub = aod::bkgds::candidateId;
  PresliceOptional<aod::JTrackDstarSubs> TracksPerDstarDataSub = aod::bkgdstar::candidateId;
  PresliceOptional<aod::JTrackLcSubs> TracksPerLcDataSub = aod::bkglc::candidateId;
  PresliceOptional<aod::JTrackB0Subs> TracksPerB0DataSub = aod::bkgb0::candidateId;
  PresliceOptional<aod::JTrackBplusSubs> TracksPerBplusDataSub = aod::bkgbplus::candidateId;
  PresliceOptional<aod::JTrackXicToXiPiPiSubs> TracksPerXicToXiPiPiDataSub = aod::bkgxictoxipipi::candidateId;
  PresliceOptional<aod::JTrackDielectronSubs> TracksPerDielectronDataSub = aod::bkgdielectron::candidateId;
  Preslice<aod::JetParticles> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  template <typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R>
  auto selectSlicer(T const& D0Slicer, U const& DplusSlicer, V const& DsSlicer, M const& DstarSlicer, N const& LcSlicer, O const& B0Slicer, P const& BplusSlicer, Q const& XicToXiPiPiSlicer, R const& DielectronSlicer)
  {
    if constexpr (jethfutilities::isD0Table<CandidateTable>()) {
      return D0Slicer;
    } else if constexpr (jethfutilities::isDplusTable<CandidateTable>()) {
      return DplusSlicer;
    } else if constexpr (jethfutilities::isDsTable<CandidateTable>()) {
      return DsSlicer;
    } else if constexpr (jethfutilities::isDstarTable<CandidateTable>()) {
      return DstarSlicer;
    } else if constexpr (jethfutilities::isLcTable<CandidateTable>()) {
      return LcSlicer;
    } else if constexpr (jethfutilities::isB0Table<CandidateTable>()) {
      return B0Slicer;
    } else if constexpr (jethfutilities::isBplusTable<CandidateTable>()) {
      return BplusSlicer;
    } else if constexpr (jethfutilities::isXicToXiPiPiTable<CandidateTable>()) {
      return XicToXiPiPiSlicer;
    } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>()) {
      return DielectronSlicer;
    } else {
      return D0Slicer;
    }
  }

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
      std::vector<int32_t> tracks;
      std::vector<int32_t> candidates;
      std::vector<int32_t> clusters;
      for (const auto& constituent : sorted_by_pt(parentSubJet2.constituents())) {
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          tracks.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidate)) {
          candidates.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
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

  template <bool isMC, bool isSubtracted, typename T, typename U, typename V, typename M, typename N>
  void jetPairing(T const& jet, U const& tracks, V const& /*candidates*/, M const& slicer, N& pairTable)
  {
    pairJetPtVec.clear();
    pairJetEnergyVec.clear();
    pairJetThetaVec.clear();
    std::vector<std::decay_t<typename U::iterator>> tracksVec;
    std::vector<std::decay_t<typename V::iterator>> candidatesVec;
    std::vector<int32_t> tracksVecIds;
    std::vector<int32_t> candidatesVecIds;
    for (auto& constituent : jet.template tracks_as<U>()) {
      if (constituent.pt() >= pairConstituentPtMin) {
        tracksVec.push_back(constituent);
        tracksVecIds.push_back(constituent.globalIndex());
      }
    }
    for (auto& candidate : jet.template candidates_as<V>()) {
      candidatesVec.push_back(candidate);
      candidatesVecIds.push_back(candidate.globalIndex());
    }
    if (tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track1Index = 0; track1Index < tracksVec.size(); track1Index++) {
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type track2Index = track1Index + 1; track2Index < tracksVec.size(); track2Index++) {
          pairJetPtVec.push_back(tracksVec.at(track1Index).pt() * tracksVec.at(track2Index).pt());
          pairJetEnergyVec.push_back(2.0 * tracksVec.at(track1Index).energy() * tracksVec.at(track2Index).energy());
          pairJetThetaVec.push_back(jetutilities::deltaR(tracksVec.at(track1Index), tracksVec.at(track2Index)));
          pairTable(jet.globalIndex(), tracksVecIds.at(track1Index), tracksVecIds.at(track2Index), -1, -1);
        }
      }
    }
    if (candidatesVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate1Index = 0; candidate1Index < candidatesVec.size(); candidate1Index++) {
        for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate2Index = candidate1Index + 1; candidate2Index < candidatesVec.size(); candidate2Index++) {
          pairJetPtVec.push_back(candidatesVec.at(candidate1Index).pt() * candidatesVec.at(candidate2Index).pt());
          auto candidate1Energy = std::sqrt((candidatesVec.at(candidate1Index).p() * candidatesVec.at(candidate1Index).p()) + (candMass * candMass));
          auto candidate2Energy = std::sqrt((candidatesVec.at(candidate2Index).p() * candidatesVec.at(candidate2Index).p()) + (candMass * candMass));
          pairJetEnergyVec.push_back(2.0 * candidate1Energy * candidate2Energy);
          pairJetThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidate1Index), candidatesVec.at(candidate2Index)));
          pairTable(jet.globalIndex(), -1, -1, candidatesVecIds.at(candidate1Index), candidatesVecIds.at(candidate2Index));
        }
      }
    }
    if (candidatesVec.size() >= 1 && tracksVec.size() >= 1) {
      for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidateIndex = 0; candidateIndex < candidatesVec.size(); candidateIndex++) { // could just directly get the candidate and tracks here but keeping it consistent with above
        for (typename std::vector<std::decay_t<typename U::iterator>>::size_type trackIndex = 0; trackIndex < tracksVec.size(); trackIndex++) {
          pairJetPtVec.push_back(candidatesVec.at(candidateIndex).pt() * tracksVec.at(trackIndex).pt());
          auto candidateEnergy = std::sqrt((candidatesVec.at(candidateIndex).p() * candidatesVec.at(candidateIndex).p()) + (candMass * candMass));
          pairJetEnergyVec.push_back(2.0 * candidateEnergy * tracksVec.at(trackIndex).energy());
          pairJetThetaVec.push_back(jetutilities::deltaR(candidatesVec.at(candidateIndex), tracksVec.at(trackIndex)));
          pairTable(jet.globalIndex(), tracksVecIds.at(trackIndex), -1, candidatesVecIds.at(candidateIndex), -1);
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
    int32_t slicerId = -1;
    if constexpr (isSubtracted) {
      auto const& candidate = jet.template candidates_first_as<V>();
      slicerId = candidate.globalIndex();
    } else {
      if constexpr (!isMC) {
        slicerId = jet.collisionId();
      } else {
        slicerId = jet.mcCollisionId();
      }
    }
    auto tracksPerCollision = tracks.sliceBy(slicer, slicerId);

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
      bool isdaughterTrack = false;
      for (auto& candidate : jet.template candidates_as<V>()) {
        if (jetcandidateutilities::isDaughterTrack(track, candidate)) {
          isdaughterTrack = true;
          break;
        }
      }
      if (isdaughterTrack) {
        continue;
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

      if (candidatesVec.size() >= 1 && tracksPerpCone1Vec.size() >= 1) {
        for (typename std::vector<std::decay_t<typename V::iterator>>::size_type candidate1Index = 0; candidate1Index < candidatesVec.size(); candidate1Index++) {
          for (typename std::vector<typename U::iterator>::size_type track2Index = 0; track2Index < tracksPerpCone1Vec.size(); track2Index++) {
            pairJetPerpCone1PtVec.push_back(candidatesVec.at(candidate1Index).pt() * tracksPerpCone1Vec.at(track2Index).pt());
            auto candidate1Energy = std::sqrt((candidatesVec.at(candidate1Index).p() * candidatesVec.at(candidate1Index).p()) + (candMass * candMass));
            pairJetPerpCone1EnergyVec.push_back(2.0 * candidate1Energy * tracksPerpCone1Vec.at(track2Index).energy());
            float dPhi = RecoDecay::constrainAngle(candidatesVec.at(candidate1Index).phi() - (tracksPerpCone1Vec.at(track2Index).phi() - (M_PI / 2.)), -M_PI);
            float dEta = candidatesVec.at(candidate1Index).eta() - tracksPerpCone1Vec.at(track2Index).eta();
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

  template <bool isSubtracted, typename T, typename U, typename V, typename M, typename N, typename O, typename P>
  void analyseCharged(T const& jet, U const& tracks, V const& candidates, M const& trackSlicer, N& outputTable, O& splittingTable, P& pairTable)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    for (auto& jetHFCandidate : jet.template candidates_as<V>()) { // should only be one at the moment
      fastjetutilities::fillTracks(jetHFCandidate, jetConstituents, jetHFCandidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidate), candMass);
    }
    nSub = jetsubstructureutilities::getNSubjettiness(jet, tracks, tracks, candidates, 2, fastjet::contrib::CA_Axes(), true, zCut, beta);
    jetReclustering<false, isSubtracted>(jet, splittingTable);
    jetPairing<false, isSubtracted>(jet, tracks, candidates, trackSlicer, pairTable);
    jetSubstructureSimple(jet, tracks, candidates);
    outputTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }

  void processChargedJetsData(typename JetTableData::iterator const& jet,
                              CandidateTable const& candidates,
                              aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, TracksPerCollision, jetSubstructureDataTable, jetSplittingsDataTable, jetPairsDataTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsData, "HF jet substructure on data", false);

  void processChargedJetsDataSub(typename JetTableDataSub::iterator const& jet,
                                 CandidateTable const& candidates,
                                 TracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks, candidates, selectSlicer(TracksPerD0DataSub, TracksPerDplusDataSub, TracksPerDsDataSub, TracksPerDstarDataSub, TracksPerLcDataSub, TracksPerB0DataSub, TracksPerBplusDataSub, TracksPerXicToXiPiPiDataSub, TracksPerDielectronDataSub), jetSubstructureDataSubTable, jetSplittingsDataSubTable, jetPairsDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsDataSub, "HF jet substructure on data", false);

  void processChargedJetsMCD(typename JetTableMCD::iterator const& jet,
                             CandidateTable const& candidates,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks, candidates, TracksPerCollision, jetSubstructureMCDTable, jetSplittingsMCDTable, jetPairsMCDTable);
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
    jetReclustering<true, false>(jet, jetSplittingsMCPTable);
    jetPairing<true, false>(jet, particles, candidates, ParticlesPerMcCollision, jetPairsMCPTable);
    jetSubstructureSimple(jet, particles, candidates);
    jetSubstructureMCPTable(energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, nSub[0], nSub[1], nSub[2], pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, angularity, leadingConstituentPt, perpConeRho);
  }
  PROCESS_SWITCH(JetSubstructureHFTask, processChargedJetsMCP, "HF jet substructure on MC particle level", false);
};
