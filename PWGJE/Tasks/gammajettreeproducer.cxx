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

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/GammaJetAnalysisTree.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "TVector2.h"

#include "CommonDataFormat/InteractionRecord.h"

#include "EventFiltering/filterTables.h"

// \struct GammaJetTreeProducer
/// \brief Task to produce a tree for gamma-jet analysis, including photons (and information of isolation) and charged and full jets
/// \author Florian Jonas <florian.jonas@cern.ch>, UC Berkeley/LBNL
/// \since 02.08.2024
///
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using selectedClusters = o2::soa::Filtered<o2::aod::JClusters>;

#include "Framework/runDataProcessing.h"

struct GammaJetTreeProducer {
  // analysis tree
  // charged jets
  // photon candidates
  Produces<aod::GjChargedJets> chargedJetsTable;
  Produces<aod::GjEvents> eventsTable;
  Produces<aod::GjGammas> gammasTable;

  HistogramRegistry mHistograms{"GammaJetTreeProducerHisto"};

  // ---------------
  // Configureables
  // ---------------

  // event cuts
  Configurable<double> mVertexCut{"vertexCut", 10.0, "apply z-vertex cut with value in cm"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<std::string>
    trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackMinPt{"trackMinPt", 0.15, "minimum track pT cut"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> isoR{"isoR", 0.4, "isolation cone radius"};

  // cluster cuts
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  // Preslice<o2::aod::JClusterTracks> perClusterMatchedTracks = o2::aod::jcluster::clusterId;

  int mRunNumber = 0;
  int eventSelection = -1;
  int trackSelection = -1;

  std::unordered_map<int32_t, int32_t> collisionMapping;
  std::vector<int> triggerMaskBits;

  void init(InitContext const&)
  {
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // create histograms
    LOG(info) << "Creating histograms";

    const o2Axis ptAxis{100, 0, 100, "p_{T} (GeV/c)"};
    const o2Axis energyAxis{100, 0, 100, "E (GeV)"};
    const o2Axis m02Axis{100, 0, 3, "m02"};

    mHistograms.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("trackPt", "pT of track", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetPt", "pT of charged jet", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetpt_vs_constpt", "pT of charged jet vs pT of constituents", o2HistType::kTH2F, {ptAxis, ptAxis});
  }

  // ---------------------
  // Helper functions
  // ---------------------
  bool isTrackSelected(const auto& track)
  {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
      return false;
    }
    if (track.pt() < trackMinPt) {
      return false;
    }

    return true;
  }

  bool isEventAccepted(const auto& collision)
  {

    if (collision.posZ() > mVertexCut) {
      return false;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return false;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return false;
    }
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return false;
    }
    return true;
  }

  double ch_iso_in_cone(const auto& cluster, aod::JetTracks const& tracks, float radius = 0.4)
  {
    double iso = 0;
    for (auto track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      // make dR function live somwhere else
      float dR = jetutilities::deltaR(cluster, track);
      if (dR < radius) {
        iso += track.pt();
      }
    }
    return iso;
  }
  double ch_perp_cone_rho(const auto& cluster, aod::JetTracks const& tracks, float radius = 0.4)
  {
    double ptSumLeft = 0;
    double ptSumRight = 0;

    double cPhi = TVector2::Phi_0_2pi(cluster.phi());

    // rotate cone left by 90 degrees
    float cPhiLeft = cPhi - TMath::Pi() / 2;
    float cPhiRight = cPhi + TMath::Pi() / 2;

    // loop over tracks
    float dRLeft, dRRight;
    for (auto track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      dRLeft = jetutilities::deltaR(cluster.eta(), cPhiLeft, track.eta(), track.phi());
      dRRight = jetutilities::deltaR(cluster.eta(), cPhiRight, track.eta(), track.phi());

      if (dRLeft < radius) {
        ptSumLeft += track.pt();
      }
      if (dRRight < radius) {
        ptSumRight += track.pt();
      }
    }

    float rho = (ptSumLeft + ptSumRight) / (2 * TMath::Pi() * radius * radius);
    return rho;
  }

  // ---------------------
  // Processing functions
  // ---------------------
  void processClearMaps(aod::JetCollisions const&)
  {
    collisionMapping.clear();
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClearMaps, "process function that clears all the maps in each dataframe", true);

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  Filter clusterDefinitionSelection = (o2::aod::jcluster::definition == mClusterDefinition);
  // Process clusters
  void processClusters(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, selectedClusters const& clusters, aod::JetTracks const& tracks)
  {
    if (!isEventAccepted(collision)) {
      return;
    }

    eventsTable(collision.multiplicity(), collision.centrality(), collision.rho(), collision.eventSel(), collision.alias_raw());
    collisionMapping[collision.globalIndex()] = eventsTable.lastIndex();

    // loop over clusters
    for (auto cluster : clusters) {

      // fill histograms
      mHistograms.fill(HIST("clusterE"), cluster.energy());

      double isoraw = ch_iso_in_cone(cluster, tracks, isoR);
      double perpconerho = ch_perp_cone_rho(cluster, tracks, isoR);

      // find closest matched track
      double dEta = 0;
      double dPhi = 0;
      // double dRMin = 100;
      double p = -1;

      // auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      // for (const auto& match : tracksofcluster) {
      //   // ask the jtracks table for track with ID trackID
      //   double dR = deltaR(cluster.eta(), cluster.phi(), match.tracks_as<o2::aod::JTracks>().Eta(), match.tracks_as<o2::aod::JTracks>().Phi());
      //   if (dR < dRMin) {
      //     dRMin = dR;
      //     dEta = cluster.eta() - match.tracks_as<o2::aod::JTracks>().eta();
      //     dPhi = TVector2::Phi_0_2pi(cluster.phi()) - TVector2::Phi_0_2pi(match.tracks_as<o2::aod::JTracks>().phi());
      //     if (abs(dPhi) > M_PI) {
      //       dPhi = 2 * M_PI - abs(dPhi);
      //     }
      //     p = match.tracks_as<o2::aod::JTracks>().p();
      //   }
      // }

      // // for compression reasons make dPhi and dEta 0 if no match is found
      // if (p == -1) {
      //   dPhi = 0;
      //   dEta = 0;
      // }

      gammasTable(eventsTable.lastIndex(), cluster.energy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(), cluster.nlm(), isoraw, perpconerho, dPhi, dEta, p);
    }

    // dummy loop over tracks
    for (auto track : tracks) {
      mHistograms.fill(HIST("trackPt"), track.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClusters, "Process EMCal clusters", true);

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  // Process charged jets
  void processChargedJets(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedJets, aod::JetTracks const&)
  {
    // event selection
    if (!isEventAccepted(collision)) {
      return;
    }

    // loop over charged jets
    for (auto jet : chargedJets) {
      if (jet.pt() < jetPtMin)
        continue;
      ushort nconst = 0;
      // loop over constituents
      for (auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
        mHistograms.fill(HIST("chjetpt_vs_constpt"), jet.pt(), constituent.pt());
        nconst++;
      }
      int32_t storedColIndex = -1;
      if (auto foundCol = collisionMapping.find(collision.globalIndex()); foundCol != collisionMapping.end()) {
        storedColIndex = foundCol->second;
      }
      chargedJetsTable(storedColIndex, jet.pt(), jet.eta(), jet.phi(), jet.energy(), jet.mass(), jet.area(), nconst);
      // fill histograms
      mHistograms.fill(HIST("chjetPt"), jet.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processChargedJets, "Process charged jets", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<GammaJetTreeProducer>(cfgc, TaskName{"gamma-jet-tree-producer"})};
  return workflow;
}
