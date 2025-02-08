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

/// \file jetderiveddatawriter.cxx
/// \brief Task to skim jet framework tables (aod::JetCollisions, aod::JetTracks, aod::JetClusters, ...)
/// while adjusting indices accordingly
///
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <MathUtils/Utils.h>
#include <algorithm>
#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  struct : ConfigurableGroup {
    Configurable<bool> performTrackSelection{"performTrackSelection", true, "only save tracks that pass one of the track selections"};
    Configurable<float> trackPtSelectionMin{"trackPtSelectionMin", 0.15, "only save tracks that have a pT larger than this pT"};
    Configurable<float> trackEtaSelectionMax{"trackEtaSelectionMax", 0.9, "only save tracks that have an eta smaller than this eta"};

    Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  } config;

  struct : ProducesGroup {
    Produces<aod::StoredJDummys> storedJDummysTable;
    Produces<aod::StoredJBCs> storedJBCsTable;
    Produces<aod::StoredJBCPIs> storedJBCParentIndexTable;
    Produces<aod::StoredJCollisions> storedJCollisionsTable;
    Produces<aod::StoredJCollisionMcInfos> storedJCollisionMcInfosTable;
    Produces<aod::StoredJCollisionPIs> storedJCollisionsParentIndexTable;
    Produces<aod::StoredJCollisionBCs> storedJCollisionsBunchCrossingIndexTable;
    Produces<aod::StoredJEMCCollisionLbs> storedJCollisionsEMCalLabelTable;
    Produces<aod::StoredJMcCollisionLbs> storedJMcCollisionsLabelTable;
    Produces<aod::StoredJMcCollisions> storedJMcCollisionsTable;
    Produces<aod::StoredJMcCollisionPIs> storedJMcCollisionsParentIndexTable;
    Produces<aod::StoredJTracks> storedJTracksTable;
    Produces<aod::StoredJTrackExtras> storedJTracksExtraTable;
    Produces<aod::StoredJEMCTracks> storedJTracksEMCalTable;
    Produces<aod::StoredJTrackPIs> storedJTracksParentIndexTable;
    Produces<aod::StoredJMcTrackLbs> storedJMcTracksLabelTable;
    Produces<aod::StoredJMcParticles> storedJMcParticlesTable;
    Produces<aod::StoredJMcParticlePIs> storedJParticlesParentIndexTable;
    Produces<aod::StoredJClusters> storedJClustersTable;
    Produces<aod::StoredJClustersCorrectedEnergies> storedJClustersCorrectedEnergiesTable;
    Produces<aod::StoredJClusterPIs> storedJClustersParentIndexTable;
    Produces<aod::StoredJClusterTracks> storedJClustersMatchedTracksTable;
    Produces<aod::StoredJMcClusterLbs> storedJMcClustersLabelTable;

    Produces<aod::StoredHfD0CollBases> storedD0CollisionsTable;
    Produces<aod::StoredJD0CollisionIds> storedD0CollisionIdsTable;
    Produces<aod::StoredHfD0Bases> storedD0sTable;
    Produces<aod::StoredHfD0Pars> storedD0ParsTable;
    Produces<aod::StoredHfD0ParEs> storedD0ParExtrasTable;
    Produces<aod::JDumD0ParDaus> storedD0ParDaughtersDummyTable;
    Produces<aod::StoredHfD0Sels> storedD0SelsTable;
    Produces<aod::StoredHfD0Mls> storedD0MlsTable;
    Produces<aod::JDumD0MlDaus> storedD0MlDughtersDummyTable;
    Produces<aod::StoredHfD0Mcs> storedD0McsTable;
    Produces<aod::StoredJD0Ids> storedD0IdsTable;
    Produces<aod::StoredHfD0McCollBases> storedD0McCollisionsTable;
    Produces<aod::StoredJD0McCollisionIds> storedD0McCollisionIdsTable;
    Produces<aod::StoredHfD0McRCollIds> storedD0McCollisionsMatchingTable;
    Produces<aod::StoredHfD0PBases> storedD0ParticlesTable;
    Produces<aod::StoredJD0PIds> storedD0ParticleIdsTable;

    Produces<aod::StoredHfDplusCollBases> storedDplusCollisionsTable;
    Produces<aod::StoredJDplusCollisionIds> storedDplusCollisionIdsTable;
    Produces<aod::StoredHfDplusBases> storedDplussTable;
    Produces<aod::StoredHfDplusPars> storedDplusParsTable;
    Produces<aod::StoredHfDplusParEs> storedDplusParExtrasTable;
    Produces<aod::JDumDplusParDaus> storedDplusParDaughtersDummyTable;
    Produces<aod::StoredHfDplusSels> storedDplusSelsTable;
    Produces<aod::StoredHfDplusMls> storedDplusMlsTable;
    Produces<aod::JDumDplusMlDaus> storedDplusMlDughtersDummyTable;
    Produces<aod::StoredHfDplusMcs> storedDplusMcsTable;
    Produces<aod::StoredJDplusIds> storedDplusIdsTable;
    Produces<aod::StoredHfDplusMcCollBases> storedDplusMcCollisionsTable;
    Produces<aod::StoredJDplusMcCollisionIds> storedDplusMcCollisionIdsTable;
    Produces<aod::StoredHfDplusMcRCollIds> storedDplusMcCollisionsMatchingTable;
    Produces<aod::StoredHfDplusPBases> storedDplusParticlesTable;
    Produces<aod::StoredJDplusPIds> storedDplusParticleIdsTable;

    Produces<aod::StoredHfLcCollBases> storedLcCollisionsTable;
    Produces<aod::StoredJLcCollisionIds> storedLcCollisionIdsTable;
    Produces<aod::StoredHfLcBases> storedLcsTable;
    Produces<aod::StoredHfLcPars> storedLcParsTable;
    Produces<aod::StoredHfLcParEs> storedLcParExtrasTable;
    Produces<aod::JDumLcParDaus> storedLcParDaughtersDummyTable;
    Produces<aod::StoredHfLcSels> storedLcSelsTable;
    Produces<aod::StoredHfLcMls> storedLcMlsTable;
    Produces<aod::JDumLcMlDaus> storedLcMlDughtersDummyTable;
    Produces<aod::StoredHfLcMcs> storedLcMcsTable;
    Produces<aod::StoredJLcIds> storedLcIdsTable;
    Produces<aod::StoredHfLcMcCollBases> storedLcMcCollisionsTable;
    Produces<aod::StoredJLcMcCollisionIds> storedLcMcCollisionIdsTable;
    Produces<aod::StoredHfLcMcRCollIds> storedLcMcCollisionsMatchingTable;
    Produces<aod::StoredHfLcPBases> storedLcParticlesTable;
    Produces<aod::StoredJLcPIds> storedLcParticleIdsTable;

    Produces<aod::StoredHfBplusCollBases> storedBplusCollisionsTable;
    Produces<aod::StoredJBplusCollisionIds> storedBplusCollisionIdsTable;
    Produces<aod::StoredHfBplusBases> storedBplussTable;
    Produces<aod::StoredHfBplusPars> storedBplusParsTable;
    Produces<aod::StoredHfBplusParEs> storedBplusParExtrasTable;
    Produces<aod::StoredHfBplusParD0s> storedBplusParD0sTable;
    Produces<aod::StoredHfBplusSels> storedBplusSelsTable;
    Produces<aod::StoredHfBplusMls> storedBplusMlsTable;
    Produces<aod::StoredHfBplusMlD0s> storedBplusMlD0sTable;
    Produces<aod::StoredHfBplusMcs> storedBplusMcsTable;
    Produces<aod::StoredJBplusIds> storedBplusIdsTable;
    Produces<aod::StoredHfBplusMcCollBases> storedBplusMcCollisionsTable;
    Produces<aod::StoredJBplusMcCollisionIds> storedBplusMcCollisionIdsTable;
    Produces<aod::StoredHfBplusMcRCollIds> storedBplusMcCollisionsMatchingTable;
    Produces<aod::StoredHfBplusPBases> storedBplusParticlesTable;
    Produces<aod::StoredJBplusPIds> storedBplusParticleIdsTable;

    Produces<aod::StoredReducedEvents> storedDielectronCollisionsTable;
    Produces<aod::StoredJDielectronCollisionIds> storedDielectronCollisionIdsTable;
    Produces<aod::StoredDielectrons> storedDielectronsTable;
    Produces<aod::StoredJDielectronIds> storedDielectronIdsTable;
    Produces<aod::StoredJDielectronMcCollisions> storedDielectronMcCollisionsTable;
    Produces<aod::StoredJDielectronMcCollisionIds> storedDielectronMcCollisionIdsTable;
    // Produces<aod::StoredHfD0McRCollIds> storedD0McCollisionsMatchingTable; //this doesnt exist for Dileptons yet
    Produces<aod::StoredJDielectronMcs> storedDielectronParticlesTable;
    Produces<aod::StoredJDielectronMcIds> storedDielectronParticleIdsTable;
  } products;

  Preslice<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>> TracksPerCollisionData = aod::jtrack::collisionId;

  Preslice<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs, aod::JMcTrackLbs>> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<aod::McCollisionsD0> D0McCollisionsPerMcCollision = aod::jd0indices::mcCollisionId;
  Preslice<aod::McCollisionsDplus> DplusMcCollisionsPerMcCollision = aod::jdplusindices::mcCollisionId;
  Preslice<aod::McCollisionsLc> LcMcCollisionsPerMcCollision = aod::jlcindices::mcCollisionId;
  Preslice<aod::McCollisionsBplus> BplusMcCollisionsPerMcCollision = aod::jbplusindices::mcCollisionId;
  Preslice<aod::McCollisionsDielectron> DielectronMcCollisionsPerMcCollision = aod::jdielectronindices::mcCollisionId;
  Preslice<aod::CandidatesD0MCP> D0ParticlesPerMcCollision = aod::jd0indices::mcCollisionId;
  Preslice<aod::CandidatesDplusMCP> DplusParticlesPerMcCollision = aod::jdplusindices::mcCollisionId;
  Preslice<aod::CandidatesLcMCP> LcParticlesPerMcCollision = aod::jlcindices::mcCollisionId;
  Preslice<aod::CandidatesBplusMCP> BplusParticlesPerMcCollision = aod::jbplusindices::mcCollisionId;
  PresliceUnsorted<aod::JEMCTracks> EMCTrackPerTrack = aod::jemctrack::trackId;

  uint32_t precisionPositionMask;
  uint32_t precisionMomentumMask;

  void init(InitContext&)
  {
    precisionPositionMask = 0xFFFFFC00; // 13 bits
    precisionMomentumMask = 0xFFFFFC00; // 13 bits  this is currently keept at 13 bits wihich gives roughly a resolution of 1/8000. This can be increased to 15 bits if really needed
  }

  template <typename T>
  bool trackSelection(T const& track)
  {
    if (config.performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept. This includes also the track pT selction.
      return false;
    }
    if (track.pt() < config.trackPtSelectionMin || std::abs(track.eta()) > config.trackEtaSelectionMax) {
      return false;
    }
    return true;
  }

  template <bool isMc, typename T>
  void storeD0(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsD0 const& D0Collisions, T const& D0s)
  {

    if (collision.isCollisionSelected()) {
      for (const auto& D0Collision : D0Collisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(D0Collision, products.storedD0CollisionsTable);
        products.storedD0CollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& D0 : D0s) {
        jethfutilities::fillHFCandidateTable<isMc>(D0, products.storedD0CollisionsTable.lastIndex(), products.storedD0sTable, products.storedD0ParsTable, products.storedD0ParExtrasTable, products.storedD0ParDaughtersDummyTable, products.storedD0SelsTable, products.storedD0MlsTable, products.storedD0MlDughtersDummyTable, products.storedD0McsTable);
        products.storedD0IdsTable(collisionMapping[collision.globalIndex()], trackMapping[D0.prong0Id()], trackMapping[D0.prong1Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeDplus(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDplus const& DplusCollisions, T const& Dpluss)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& DplusCollision : DplusCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(DplusCollision, products.storedDplusCollisionsTable);
        products.storedDplusCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& Dplus : Dpluss) {
        jethfutilities::fillHFCandidateTable<isMc>(Dplus, products.storedDplusCollisionsTable.lastIndex(), products.storedDplussTable, products.storedDplusParsTable, products.storedDplusParExtrasTable, products.storedDplusParDaughtersDummyTable, products.storedDplusSelsTable, products.storedDplusMlsTable, products.storedDplusMlDughtersDummyTable, products.storedDplusMcsTable);
        products.storedDplusIdsTable(collisionMapping[collision.globalIndex()], trackMapping[Dplus.prong0Id()], trackMapping[Dplus.prong1Id()], trackMapping[Dplus.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeLc(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsLc const& LcCollisions, T const& Lcs)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& LcCollision : LcCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(LcCollision, products.storedLcCollisionsTable);
        products.storedLcCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& Lc : Lcs) {
        jethfutilities::fillHFCandidateTable<isMc>(Lc, products.storedLcCollisionsTable.lastIndex(), products.storedLcsTable, products.storedLcParsTable, products.storedLcParExtrasTable, products.storedLcParDaughtersDummyTable, products.storedLcSelsTable, products.storedLcMlsTable, products.storedLcMlDughtersDummyTable, products.storedLcMcsTable);
        products.storedLcIdsTable(collisionMapping[collision.globalIndex()], trackMapping[Lc.prong0Id()], trackMapping[Lc.prong1Id()], trackMapping[Lc.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeBplus(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsBplus const& BplusCollisions, T const& Bpluss)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& BplusCollision : BplusCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(BplusCollision, products.storedBplusCollisionsTable);
        products.storedBplusCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& Bplus : Bpluss) {
        jethfutilities::fillHFCandidateTable<isMc>(Bplus, products.storedBplusCollisionsTable.lastIndex(), products.storedBplussTable, products.storedBplusParsTable, products.storedBplusParExtrasTable, products.storedBplusParD0sTable, products.storedBplusSelsTable, products.storedBplusMlsTable, products.storedBplusMlD0sTable, products.storedBplusMcsTable);
        products.storedBplusIdsTable(collisionMapping[collision.globalIndex()], trackMapping[Bplus.prong0Id()], trackMapping[Bplus.prong1Id()], trackMapping[Bplus.prong2Id()]);
      }
    }
  }

  void processDummyTable(aod::JDummys const&)
  {
    products.storedJDummysTable(1);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDummyTable, "write out dummy output table", true);

  std::vector<int32_t> collisionMapping;
  std::vector<int32_t> bcMapping;
  std::vector<int32_t> trackMapping;
  std::vector<int32_t> mcCollisionMapping;
  std::vector<int32_t> particleMapping;
  std::vector<int32_t> d0McCollisionMapping;
  std::vector<int32_t> dplusMcCollisionMapping;
  std::vector<int32_t> lcMcCollisionMapping;
  std::vector<int32_t> bplusMcCollisionMapping;
  std::vector<int32_t> dielectronMcCollisionMapping;

  void processBCs(soa::Join<aod::JCollisions, aod::JCollisionBCs, aod::JCollisionSelections> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const& bcs)
  {
    std::vector<int32_t> bcIndicies;
    bcMapping.clear();
    bcMapping.resize(bcs.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp(), bc.alias_raw(), bc.selection_raw());
          products.storedJBCParentIndexTable(bc.bcId());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping[bc.globalIndex()] = products.storedJBCsTable.lastIndex();
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBCs, "write out output tables for Bunch crossings", true);

  void processColllisons(soa::Join<aod::JCollisions, aod::JCollisionMcInfos, aod::JCollisionPIs, aod::JCollisionBCs, aod::JCollisionSelections> const& collisions)
  {
    collisionMapping.clear();
    collisionMapping.resize(collisions.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {

        products.storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.trackOccupancyInTimeRange(), collision.eventSel(), collision.alias_raw(), collision.triggerSel());
        collisionMapping[collision.globalIndex()] = products.storedJCollisionsTable.lastIndex();
        products.storedJCollisionMcInfosTable(collision.weight(), collision.subGeneratorId());
        products.storedJCollisionsParentIndexTable(collision.collisionId());
        if (doprocessBCs) {
          products.storedJCollisionsBunchCrossingIndexTable(bcMapping[collision.bcId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processColllisons, "write out output tables for collisions", true);

  void processTracks(soa::Join<aod::JCollisions, aod::JCollisionSelections> const& collisions, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& tracks)
  {
    trackMapping.clear();
    trackMapping.resize(tracks.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        const auto tracksPerCollision = tracks.sliceBy(TracksPerCollisionData, collision.globalIndex());
        for (const auto& track : tracksPerCollision) {
          if (!trackSelection(track)) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept. This includes also the track pT selction.
            continue;
          }
          products.storedJTracksTable(collisionMapping[collision.globalIndex()], o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
          products.storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaX(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigma1Pt(), precisionMomentumMask));
          products.storedJTracksParentIndexTable(track.trackId());
          trackMapping[track.globalIndex()] = products.storedJTracksTable.lastIndex();
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processTracks, "write out output tables for tracks", true);

  void processClusters(soa::Join<aod::JCollisions, aod::JEMCCollisionLbs, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::JEMCTracks const& emcTracks, soa::Join<aod::JClusters, aod::JClustersCorrectedEnergies, aod::JClusterPIs, aod::JClusterTracks> const& clusters)
  {
    if (collision.isCollisionSelected()) {
      products.storedJCollisionsEMCalLabelTable(collision.isAmbiguous(), collision.isEmcalReadout());
      for (const auto& cluster : clusters) {
        products.storedJClustersTable(collisionMapping[collision.globalIndex()], cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                      cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                      cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
        products.storedJClustersCorrectedEnergiesTable(cluster.energyCorrectedOneTrack1(), cluster.energyCorrectedOneTrack2(), cluster.energyCorrectedAllTracks1(), cluster.energyCorrectedAllTracks2());
        products.storedJClustersParentIndexTable(cluster.clusterId());

        std::vector<int32_t> clusterStoredJTrackIDs;
        for (const auto& clusterTrack : cluster.matchedTracks_as<aod::JTracks>()) {
          clusterStoredJTrackIDs.push_back(trackMapping[clusterTrack.globalIndex()]);
          auto emcTracksPerTrack = emcTracks.sliceBy(EMCTrackPerTrack, clusterTrack.globalIndex());
          auto emcTrackPerTrack = emcTracksPerTrack.iteratorAt(0);
          products.storedJTracksEMCalTable(trackMapping[clusterTrack.globalIndex()], emcTrackPerTrack.etaEmcal(), emcTrackPerTrack.phiEmcal());
        }
        products.storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processClusters, "write out output tables for clusters", false);

  //!!!!!!!!!! need to add the hadronic corrected energy and delete the new dummy process function

  void processD0Data(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0Data const& D0s)
  {
    storeD0<false>(collision, tracks, D0Collisions, D0s);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0Data, "write out data output tables for D0", false);

  void processD0MCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0MCD const& D0s)
  {
    storeD0<true>(collision, tracks, D0Collisions, D0s);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0MCD, "write out mcd output tables for D0", false);

  void processDplusData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDplus const& DplusCollisions, aod::CandidatesDplusData const& Dpluss)
  {
    storeDplus<false>(collision, tracks, DplusCollisions, Dpluss);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusData, "write out data output tables for Dplus", false);

  void processDplusMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDplus const& DplusCollisions, aod::CandidatesDplusMCD const& Dpluss)
  {
    storeDplus<true>(collision, tracks, DplusCollisions, Dpluss);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMCD, "write out mcd output tables for Dplus", false);

  void processLcData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcData const& Lcs)
  {
    storeLc<false>(collision, tracks, LcCollisions, Lcs);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcData, "write out data output tables for Lc", false);

  void processLcMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcMCD const& Lcs)
  {
    storeLc<true>(collision, tracks, LcCollisions, Lcs);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcMCD, "write out mcd output tables for Lc", false);

  void processBplusData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsBplus const& BplusCollisions, aod::CandidatesBplusData const& Bpluss)
  {
    storeBplus<false>(collision, tracks, BplusCollisions, Bpluss);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusData, "write out data output tables for bplus", false);

  void processBplusMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsBplus const& BplusCollisions, aod::CandidatesBplusMCD const& Bpluss)
  {
    storeBplus<true>(collision, tracks, BplusCollisions, Bpluss);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMCD, "write out mcd output tables for bplus", false);

  void processDielectron(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDielectron const& DielectronCollisions, aod::CandidatesDielectronData const& Dielectrons)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& DielectronCollision : DielectronCollisions) { // should only ever be one
        jetdqutilities::fillDielectronCollisionTable(DielectronCollision, products.storedDielectronCollisionsTable);
        products.storedDielectronCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& Dielectron : Dielectrons) {
        jetdqutilities::fillDielectronCandidateTable(Dielectron, products.storedDielectronCollisionsTable.lastIndex(), products.storedDielectronsTable);
        products.storedDielectronIdsTable(collisionMapping[collision.globalIndex()], trackMapping[Dielectron.prong0Id()], trackMapping[Dielectron.prong1Id()]);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDielectron, "write out data output tables for dielectrons", false);

  void processMcCollisions(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs, aod::JMcCollisionSelections> const& mcCollisions)
  {
    mcCollisionMapping.clear();
    mcCollisionMapping.resize(mcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {

        products.storedJMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.weight(), mcCollision.subGeneratorId());
        products.storedJMcCollisionsParentIndexTable(mcCollision.mcCollisionId());
        mcCollisionMapping[mcCollision.globalIndex()] = products.storedJMcCollisionsTable.lastIndex();
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processMcCollisions, "write out mcCollision output tables", false);

  void processMcParticles(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles)
  {
    particleMapping.clear();
    particleMapping.resize(particles.size(), -1);
    int particleTableIndex = 0;
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {

        const auto particlesPerMcCollision = particles.sliceBy(ParticlesPerMcCollision, mcCollision.globalIndex());

        for (auto particle : particlesPerMcCollision) {
          particleMapping[particle.globalIndex()] = particleTableIndex;
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {

          std::vector<int32_t> mothersIds;
          if (particle.has_mothers()) {
            auto mothersIdTemps = particle.mothersIds();
            for (auto mothersIdTemp : mothersIdTemps) {
              mothersIds.push_back(particleMapping[mothersIdTemp]);
            }
          }
          int daughtersIds[2] = {-1, -1};
          auto i = 0;
          if (particle.has_daughters()) {
            for (auto daughterId : particle.daughtersIds()) {
              if (i > 1) {
                break;
              }
              daughtersIds[i] = particleMapping[daughterId];
              i++;
            }
          }
          products.storedJMcParticlesTable(mcCollisionMapping[mcCollision.globalIndex()], o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersIds, daughtersIds);
          products.storedJParticlesParentIndexTable(particle.mcParticleId());
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processMcParticles, "write out mcParticle output tables", false);

  void processD0MCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsD0 const& D0McCollisions, aod::CandidatesD0MCP const& D0Particles)
  {
    d0McCollisionMapping.clear();
    d0McCollisionMapping.resize(D0McCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto d0McCollisionsPerMcCollision = D0McCollisions.sliceBy(D0McCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& d0McCollisionPerMcCollision : d0McCollisionsPerMcCollision) {
          jethfutilities::fillHFMcCollisionTable(d0McCollisionPerMcCollision, products.storedD0McCollisionsTable);
          products.storedD0McCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          d0McCollisionMapping[d0McCollisionPerMcCollision.globalIndex()] = products.storedD0McCollisionsTable.lastIndex();
        }
        const auto d0ParticlesPerMcCollision = D0Particles.sliceBy(D0ParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& D0Particle : d0ParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(D0Particle, products.storedD0McCollisionsTable.lastIndex(), products.storedD0ParticlesTable);
          products.storedD0ParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[D0Particle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0MCP, "write out D0 mcp output tables", false);

  void processDplusMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsDplus const& DplusMcCollisions, aod::CandidatesDplusMCP const& DplusParticles)
  {
    dplusMcCollisionMapping.clear();
    dplusMcCollisionMapping.resize(DplusMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto dplusMcCollisionsPerMcCollision = DplusMcCollisions.sliceBy(DplusMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& dplusMcCollisionPerMcCollision : dplusMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(dplusMcCollisionPerMcCollision, products.storedDplusMcCollisionsTable);
          products.storedDplusMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          dplusMcCollisionMapping[dplusMcCollisionPerMcCollision.globalIndex()] = products.storedDplusMcCollisionsTable.lastIndex();
        }
        const auto dplusParticlesPerMcCollision = DplusParticles.sliceBy(DplusParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& DplusParticle : dplusParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(DplusParticle, products.storedDplusMcCollisionsTable.lastIndex(), products.storedDplusParticlesTable);
          products.storedDplusParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DplusParticle.mcParticleId()]);
        }
      }
    }
  }

  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMCP, "write out Dplus mcp output tables", false);

  void processLcMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsLc const& LcMcCollisions, aod::CandidatesLcMCP const& LcParticles)
  {
    lcMcCollisionMapping.clear();
    lcMcCollisionMapping.resize(LcMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto lcMcCollisionsPerMcCollision = LcMcCollisions.sliceBy(LcMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& lcMcCollisionPerMcCollision : lcMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(lcMcCollisionPerMcCollision, products.storedLcMcCollisionsTable);
          products.storedLcMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          lcMcCollisionMapping[lcMcCollisionPerMcCollision.globalIndex()] = products.storedLcMcCollisionsTable.lastIndex();
        }
        const auto lcParticlesPerMcCollision = LcParticles.sliceBy(LcParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& LcParticle : lcParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(LcParticle, products.storedLcMcCollisionsTable.lastIndex(), products.storedLcParticlesTable);
          products.storedLcParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[LcParticle.mcParticleId()]);
        }
      }
    }
  }

  PROCESS_SWITCH(JetDerivedDataWriter, processLcMCP, "write out Lc mcp output tables", false);

  void processBplusMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsBplus const& BplusMcCollisions, aod::CandidatesBplusMCP const& BplusParticles)
  {
    bplusMcCollisionMapping.clear();
    bplusMcCollisionMapping.resize(BplusMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto bplusMcCollisionsPerMcCollision = BplusMcCollisions.sliceBy(BplusMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& bplusMcCollisionPerMcCollision : bplusMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(bplusMcCollisionPerMcCollision, products.storedBplusMcCollisionsTable);
          products.storedBplusMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          bplusMcCollisionMapping[bplusMcCollisionPerMcCollision.globalIndex()] = products.storedBplusMcCollisionsTable.lastIndex();
        }
        const auto bplusParticlesPerMcCollision = BplusParticles.sliceBy(BplusParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& BplusParticle : bplusParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(BplusParticle, products.storedBplusMcCollisionsTable.lastIndex(), products.storedBplusParticlesTable);
          products.storedBplusParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[BplusParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMCP, "write out Bplus mcp output tables", false);

  void processDielectronMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, aod::JMcParticles const&, aod::McCollisionsDielectron const& DielectronMcCollisions, aod::CandidatesDielectronMCP const& DielectronParticles)
  {
    if (mcCollision.isMcCollisionSelected()) {

      const auto dielectronMcCollisionsPerMcCollision = DielectronMcCollisions.sliceBy(DielectronMcCollisionsPerMcCollision, mcCollision.globalIndex());
      for (const auto& dielectronMcCollisionPerMcCollision : dielectronMcCollisionsPerMcCollision) { // should only ever be one
        jetdqutilities::fillDielectronMcCollisionTable(dielectronMcCollisionPerMcCollision, products.storedDielectronMcCollisionsTable);
        products.storedDielectronMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
      }
      for (const auto& DielectronParticle : DielectronParticles) {
        jetdqutilities::fillDielectronCandidateMcTable(DielectronParticle, products.storedDielectronMcCollisionsTable.lastIndex(), products.storedDielectronParticlesTable);
        std::vector<int32_t> DielectronMothersIds;
        int DielectronDaughtersId[2];
        if (DielectronParticle.has_mothers()) {
          for (auto const& DielectronMother : DielectronParticle.template mothers_as<aod::JMcParticles>()) {
            DielectronMothersIds.push_back(particleMapping[DielectronMother.globalIndex()]);
          }
        }
        auto i = 0;
        if (DielectronParticle.has_daughters()) {
          for (auto const& DielectronDaughter : DielectronParticle.template daughters_as<aod::JMcParticles>()) {
            if (i > 1) {
              break;
            }
            DielectronDaughtersId[i] = particleMapping[DielectronDaughter.globalIndex()];
            i++;
          }
        }
        products.storedDielectronParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DielectronParticle.mcParticleId()], DielectronMothersIds, DielectronDaughtersId);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDielectronMCP, "write out Dielectron mcp output tables", false);

  void processColllisonsMcCollisionLabel(soa::Join<aod::JCollisions, aod::JMcCollisionLbs, aod::JCollisionSelections>::iterator const& collision)
  {
    if (collision.isCollisionSelected()) {
      products.storedJMcCollisionsLabelTable(mcCollisionMapping[collision.mcCollisionId()]);
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processColllisonsMcCollisionLabel, "write out collision mcCollision label output tables", false);

  void processTracksMcParticleLabel(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, soa::Join<aod::JTracks, aod::JMcTrackLbs> const& tracks)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& track : tracks) {
        if (!trackSelection(track)) {
          continue;
        }
        if (track.has_mcParticle()) {
          products.storedJMcTracksLabelTable(particleMapping[track.mcParticleId()]);
        } else {
          products.storedJMcTracksLabelTable(-1);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processTracksMcParticleLabel, "write out track mcParticle label output tables", false);

  void processClusterMcLabel(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, soa::Join<aod::JClusters, aod::JMcClusterLbs> const& clusters)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& cluster : clusters) {
        std::vector<int32_t> clusterStoredJParticleIDs;
        for (const auto& clusterParticleId : cluster.mcParticlesIds()) {
          clusterStoredJParticleIDs.push_back(particleMapping[clusterParticleId]);
        }
        std::vector<float> amplitudeA;
        auto amplitudeASpan = cluster.amplitudeA();
        std::copy(amplitudeASpan.begin(), amplitudeASpan.end(), std::back_inserter(amplitudeA));
        products.storedJMcClustersLabelTable(clusterStoredJParticleIDs, amplitudeA);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processClusterMcLabel, "write out cluster mc label output tables", false);

  void processD0McCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsD0, aod::HfD0McRCollIds> const& D0McCollisions, aod::CollisionsD0 const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& D0McCollision : D0McCollisions) { // should just be one
        std::vector<int32_t> d0CollisionIDs;
        for (auto const& d0CollisionPerMcCollision : D0McCollision.hfCollBases_as<aod::CollisionsD0>()) {
          d0CollisionIDs.push_back(d0McCollisionMapping[d0CollisionPerMcCollision.globalIndex()]);
        }
        products.storedD0McCollisionsMatchingTable(d0CollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0McCollisionMatch, "write out D0 McCollision collision label output tables", false);

  void processDplusMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsDplus, aod::HfDplusMcRCollIds> const& DplusMcCollisions, aod::CollisionsDplus const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& DplusMcCollision : DplusMcCollisions) { // should just be one
        std::vector<int32_t> dplusCollisionIDs;
        for (auto const& dplusCollisionPerMcCollision : DplusMcCollision.hfCollBases_as<aod::CollisionsDplus>()) {
          dplusCollisionIDs.push_back(dplusMcCollisionMapping[dplusCollisionPerMcCollision.globalIndex()]);
        }
        products.storedDplusMcCollisionsMatchingTable(dplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMcCollisionMatch, "write out Dplus McCollision collision label output tables", false);

  void processLcMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsLc, aod::HfLcMcRCollIds> const& LcMcCollisions, aod::CollisionsLc const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& LcMcCollision : LcMcCollisions) { // should just be one
        std::vector<int32_t> lcCollisionIDs;
        for (auto const& lcCollisionPerMcCollision : LcMcCollision.hfCollBases_as<aod::CollisionsLc>()) {
          lcCollisionIDs.push_back(lcMcCollisionMapping[lcCollisionPerMcCollision.globalIndex()]);
        }
        products.storedLcMcCollisionsMatchingTable(lcCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcMcCollisionMatch, "write out Lc McCollision collision label output tables", false);

  void processBplusMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsBplus, aod::HfBplusMcRCollIds> const& BplusMcCollisions, aod::CollisionsBplus const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& BplusMcCollision : BplusMcCollisions) { // should just be one
        std::vector<int32_t> bplusCollisionIDs;
        for (auto const& bplusCollisionPerMcCollision : BplusMcCollision.hfCollBases_as<aod::CollisionsBplus>()) {
          bplusCollisionIDs.push_back(bplusMcCollisionMapping[bplusCollisionPerMcCollision.globalIndex()]);
        }
        products.storedBplusMcCollisionsMatchingTable(bplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMcCollisionMatch, "write out Bplus McCollision collision label output tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
