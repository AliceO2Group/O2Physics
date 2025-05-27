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

// C++ system headers first
#include <string>
#include <unordered_map>
#include <vector>

// Framework and other headers after
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
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using emcClusters = o2::soa::Join<o2::aod::JClusters, o2::aod::JClusterTracks>;

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
  Configurable<float> isoR{"isoR", 0.4, "isolation cone radius"};
  Configurable<float> perpConeJetR{"perpConeJetR", 0.4, "perpendicular cone radius used to calculate perp cone rho for jet"};
  Configurable<float> trackMatchingEoverP{"trackMatchingEoverP", 2.0, "closest track is required to have E/p < value"};
  Configurable<float> minClusterETrigger{"minClusterETrigger", 0.0, "minimum cluster energy to trigger"};

  int mRunNumber = 0;
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  std::unordered_map<int32_t, int32_t> collisionMapping;
  std::vector<int> triggerMaskBits;

  void init(InitContext const&)
  {
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // create histograms
    LOG(info) << "Creating histograms";

    const o2Axis ptAxis{100, 0, 200, "p_{T} (GeV/c)"};
    const o2Axis energyAxis{100, 0, 100, "E (GeV)"};
    const o2Axis m02Axis{100, 0, 3, "m02"};
    const o2Axis etaAxis{100, -1, 1, "#eta"};
    const o2Axis phiAxis{100, 0, 2 * TMath::Pi(), "#phi"};
    const o2Axis occupancyAxis{300, 0, 30000, "occupancy"};
    mHistograms.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("trackPt", "pT of track", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetPt", "pT of charged jet", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetPtEtaPhi", "pT of charged jet", o2HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
    mHistograms.add("chjetpt_vs_constpt", "pT of charged jet vs pT of constituents", o2HistType::kTH2F, {ptAxis, ptAxis});

    // track QA THnSparse
    mHistograms.add("trackPtEtaPhi", "Track QA", o2HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
    mHistograms.add("trackPtEtaOccupancy", "Track QA vs occupancy", o2HistType::kTHnSparseF, {ptAxis, etaAxis, occupancyAxis});
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

  int getStoredColIndex(const auto& collision)
  {
    int32_t storedColIndex = -1;
    if (auto foundCol = collisionMapping.find(collision.globalIndex()); foundCol != collisionMapping.end()) {
      storedColIndex = foundCol->second;
    }
    return storedColIndex;
  }

  bool isEventAccepted(const auto& collision, const auto& clusters)
  {

    if (collision.posZ() > mVertexCut) {
      return false;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return false;
    }
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return false;
    }
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return false;
    }

    // Check if event contains a cluster with energy > minClusterETrigger
    for (auto cluster : clusters) {
      if (cluster.energy() > minClusterETrigger) {
        return true;
      }
    }
    return false;
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

  void runTrackQA(const auto& collision, aod::JetTracks const& tracks)
  {
    for (auto track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      mHistograms.fill(HIST("trackPt"), track.pt());
      mHistograms.fill(HIST("trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      mHistograms.fill(HIST("trackPtEtaOccupancy"), track.pt(), track.eta(), collision.trackOccupancyInTimeRange());
    }
  }

  double ch_perp_cone_rho(const auto& object, aod::JetTracks const& tracks, float radius = 0.4)
  {
    double ptSumLeft = 0;
    double ptSumRight = 0;

    double cPhi = TVector2::Phi_0_2pi(object.phi());

    // rotate cone left by 90 degrees
    float cPhiLeft = cPhi - TMath::Pi() / 2;
    float cPhiRight = cPhi + TMath::Pi() / 2;

    // loop over tracks
    float dRLeft, dRRight;
    for (auto track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      dRLeft = jetutilities::deltaR(object.eta(), cPhiLeft, track.eta(), track.phi());
      dRRight = jetutilities::deltaR(object.eta(), cPhiRight, track.eta(), track.phi());

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
  // WARNING: This function always has to run first in the processing chain
  void processClearMaps(aod::JetCollisions const&)
  {
    collisionMapping.clear();
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClearMaps, "process function that clears all the maps in each dataframe", true);

  // WARNING: This function always has to run second in the processing chain
  void processEvent(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, emcClusters const& clusters)
  {
    if (!isEventAccepted(collision, clusters)) {
      return;
    }

    eventsTable(collision.multiplicity(), collision.centrality(), collision.rho(), collision.eventSel(), collision.trackOccupancyInTimeRange(), collision.alias_raw());
    collisionMapping[collision.globalIndex()] = eventsTable.lastIndex();
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processEvent, "Process event", true);

  // ---------------------
  // Processing functions can be safely added below this line
  // ---------------------

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  PresliceUnsorted<aod::JEMCTracks> EMCTrackPerTrack = aod::jemctrack::trackId;
  // Process clusters
  void processClusters(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, emcClusters const& clusters, aod::JetTracks const& tracks, aod::JEMCTracks const& emctracks)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;

    // eventsTable(collision.multiplicity(), collision.centrality(), collision.rho(), collision.eventSel(), collision.trackOccupancyInTimeRange(), collision.alias_raw());
    // collisionMapping[collision.globalIndex()] = eventsTable.lastIndex();

    // loop over tracks one time for QA
    runTrackQA(collision, tracks);

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

      // do track matching
      auto tracksofcluster = cluster.matchedTracks_as<aod::JetTracks>();
      for (auto track : tracksofcluster) {
        if (!isTrackSelected(track)) {
          continue;
        }
        auto emcTracksPerTrack = emctracks.sliceBy(EMCTrackPerTrack, track.globalIndex());
        auto emcTrack = emcTracksPerTrack.iteratorAt(0);
        // find closest track that still has E/p < trackMatchingEoverP
        if (cluster.energy() / track.p() > trackMatchingEoverP) {
          continue;
        } else {
          dEta = cluster.eta() - emcTrack.etaEmcal();
          dPhi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(emcTrack.phiEmcal(), -M_PI) - RecoDecay::constrainAngle(cluster.phi(), -M_PI), -M_PI);
          p = track.p();
          break;
        }
      }
      gammasTable(storedColIndex, cluster.energy(), cluster.definition(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(), cluster.nlm(), isoraw, perpconerho, dPhi, dEta, p);
    }

    // dummy loop over tracks
    for (auto track : tracks) {
      mHistograms.fill(HIST("trackPt"), track.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClusters, "Process EMCal clusters", true);

  Filter jetCuts = aod::jet::pt > jetPtMin;
  // Process charged jets
  void processChargedJets(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedJets, aod::JetTracks const& tracks)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;
    float leadingTrackPt = 0;
    ushort nconst = 0;
    // loop over charged jets
    for (auto jet : chargedJets) {
      if (jet.pt() < jetPtMin)
        continue;
      nconst = 0;
      leadingTrackPt = 0;
      // loop over constituents
      for (auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
        mHistograms.fill(HIST("chjetpt_vs_constpt"), jet.pt(), constituent.pt());
        nconst++;
        if (constituent.pt() > leadingTrackPt) {
          leadingTrackPt = constituent.pt();
        }
      }

      // calculate perp cone rho
      double perpconerho = ch_perp_cone_rho(jet, tracks, perpConeJetR);
      mHistograms.fill(HIST("chjetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      chargedJetsTable(storedColIndex, jet.pt(), jet.eta(), jet.phi(), jet.r(), jet.energy(), jet.mass(), jet.area(), leadingTrackPt, perpconerho, nconst);
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
