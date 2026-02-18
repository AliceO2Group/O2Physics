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

#include "JetDerivedDataUtilities.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"
#include "PWGJE/DataModel/JetReducedDataSelector.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/detail/TypeTruncation.h>

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  struct : ConfigurableGroup {
    Configurable<bool> performTrackSelection{"performTrackSelection", true, "only save tracks that pass one of the track selections"};
    Configurable<float> trackPtSelectionMin{"trackPtSelectionMin", 0.15, "only save tracks that have a pT larger than this pT"};
    Configurable<float> trackEtaSelectionMax{"trackEtaSelectionMax", 0.9, "only save tracks that have an eta smaller than this eta"};
    Configurable<bool> savePartonLevelInfo{"savePartonLevelInfo", true, "save parton level info at MCP level"};

  } config;

  struct : ProducesGroup {
    Produces<aod::StoredJDummys> storedJDummysTable;
    Produces<aod::StoredJBCs> storedJBCsTable;
    Produces<aod::StoredJBCPIs> storedJBCParentIndexTable;
    Produces<aod::StoredJCollisions> storedJCollisionsTable;
    Produces<aod::StoredJCollisionUPCs> storedJCollisionUPCsTable;
    Produces<aod::StoredJCollisionMcInfos> storedJCollisionMcInfosTable;
    Produces<aod::StoredJCollisionPIs> storedJCollisionsParentIndexTable;
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
    Produces<aod::StoredJClusterPIs> storedJClustersParentIndexTable;
    Produces<aod::StoredJClusterTracks> storedJClustersMatchedTracksTable;
    Produces<aod::StoredJMcClusterLbs> storedJMcClustersLabelTable;

    struct : ProducesGroup {
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
    } productsD0;

    struct : ProducesGroup {
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
    } productsDplus;

    struct : ProducesGroup {
      Produces<aod::StoredHfDsCollBases> storedDsCollisionsTable;
      Produces<aod::StoredJDsCollisionIds> storedDsCollisionIdsTable;
      Produces<aod::StoredHfDsBases> storedDssTable;
      Produces<aod::StoredHfDsPars> storedDsParsTable;
      Produces<aod::StoredHfDsParEs> storedDsParExtrasTable;
      Produces<aod::JDumDsParDaus> storedDsParDaughtersDummyTable;
      Produces<aod::StoredHfDsSels> storedDsSelsTable;
      Produces<aod::StoredHfDsMls> storedDsMlsTable;
      Produces<aod::JDumDsMlDaus> storedDsMlDughtersDummyTable;
      Produces<aod::StoredHfDsMcs> storedDsMcsTable;
      Produces<aod::StoredJDsIds> storedDsIdsTable;
      Produces<aod::StoredHfDsMcCollBases> storedDsMcCollisionsTable;
      Produces<aod::StoredJDsMcCollisionIds> storedDsMcCollisionIdsTable;
      Produces<aod::StoredHfDsMcRCollIds> storedDsMcCollisionsMatchingTable;
      Produces<aod::StoredHfDsPBases> storedDsParticlesTable;
      Produces<aod::StoredJDsPIds> storedDsParticleIdsTable;
    } productsDs;

    struct : ProducesGroup {
      Produces<aod::StoredHfDstarCollBases> storedDstarCollisionsTable;
      Produces<aod::StoredJDstarCollisionIds> storedDstarCollisionIdsTable;
      Produces<aod::StoredHfDstarBases> storedDstarsTable;
      Produces<aod::StoredHfDstarPars> storedDstarParsTable;
      Produces<aod::JDumDstarParEs> storedDstarParExtrasDummyTable;
      Produces<aod::StoredHfDstarParD0s> storedDstarParD0sTable;
      Produces<aod::StoredHfDstarSels> storedDstarSelsTable;
      Produces<aod::StoredHfDstarMls> storedDstarMlsTable;
      Produces<aod::JDumDstarMlDaus> storedDstarMlDaughtersDummyTable;
      Produces<aod::StoredHfDstarMcs> storedDstarMcsTable;
      Produces<aod::StoredJDstarIds> storedDstarIdsTable;
      Produces<aod::StoredHfDstarMcCollBases> storedDstarMcCollisionsTable;
      Produces<aod::StoredJDstarMcCollisionIds> storedDstarMcCollisionIdsTable;
      Produces<aod::StoredHfDstarMcRCollIds> storedDstarMcCollisionsMatchingTable;
      Produces<aod::StoredHfDstarPBases> storedDstarParticlesTable;
      Produces<aod::StoredJDstarPIds> storedDstarParticleIdsTable;
    } productsDstar;

    struct : ProducesGroup {
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
    } productsLc;

    struct : ProducesGroup {
      Produces<aod::StoredHfB0CollBases> storedB0CollisionsTable;
      Produces<aod::StoredJB0CollisionIds> storedB0CollisionIdsTable;
      Produces<aod::StoredHfB0Bases> storedB0sTable;
      Produces<aod::StoredHfB0Pars> storedB0ParsTable;
      Produces<aod::StoredHfB0ParEs> storedB0ParExtrasTable;
      Produces<aod::StoredHfB0ParDpluss> storedB0ParDplussTable;
      Produces<aod::StoredHfB0Sels> storedB0SelsTable;
      Produces<aod::StoredHfB0Mls> storedB0MlsTable;
      Produces<aod::StoredHfB0MlDpluss> storedB0MlDplussTable;
      Produces<aod::StoredHfB0Mcs> storedB0McsTable;
      Produces<aod::StoredJB0Ids> storedB0IdsTable;
      Produces<aod::StoredHfB0McCollBases> storedB0McCollisionsTable;
      Produces<aod::StoredJB0McCollisionIds> storedB0McCollisionIdsTable;
      Produces<aod::StoredHfB0McRCollIds> storedB0McCollisionsMatchingTable;
      Produces<aod::StoredHfB0PBases> storedB0ParticlesTable;
      Produces<aod::StoredJB0PIds> storedB0ParticleIdsTable;
    } productsB0;

    struct : ProducesGroup {
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
    } productsBplus;

    struct : ProducesGroup {
      Produces<aod::StoredHfXicToXiPiPiCollBases> storedXicToXiPiPiCollisionsTable;
      Produces<aod::StoredJXicToXiPiPiCollisionIds> storedXicToXiPiPiCollisionIdsTable;
      Produces<aod::StoredHfXicToXiPiPiBases> storedXicToXiPiPisTable;
      Produces<aod::StoredHfXicToXiPiPiPars> storedXicToXiPiPiParsTable;
      Produces<aod::StoredHfXicToXiPiPiParEs> storedXicToXiPiPiParExtrasTable;
      Produces<aod::JDumXicToXiPiPiParDaus> storedXicToXiPiPiParDaughtersDummyTable;
      Produces<aod::StoredHfXicToXiPiPiSels> storedXicToXiPiPiSelsTable;
      Produces<aod::StoredHfXicToXiPiPiMls> storedXicToXiPiPiMlsTable;
      Produces<aod::JDumXicToXiPiPiMlDaus> storedXicToXiPiPiMlDughtersDummyTable;
      Produces<aod::StoredHfXicToXiPiPiMcs> storedXicToXiPiPiMcsTable;
      Produces<aod::StoredJXicToXiPiPiIds> storedXicToXiPiPiIdsTable;
      Produces<aod::StoredHfXicToXiPiPiMcCollBases> storedXicToXiPiPiMcCollisionsTable;
      Produces<aod::StoredJXicToXiPiPiMcCollisionIds> storedXicToXiPiPiMcCollisionIdsTable;
      Produces<aod::StoredHfXicToXiPiPiMcRCollIds> storedXicToXiPiPiMcCollisionsMatchingTable;
      Produces<aod::StoredHfXicToXiPiPiPBases> storedXicToXiPiPiParticlesTable;
      Produces<aod::StoredJXicToXiPiPiPIds> storedXicToXiPiPiParticleIdsTable;
    } productsXicToXiPiPi;

    struct : ProducesGroup {
      Produces<aod::StoredReducedEvents> storedDielectronCollisionsTable;
      Produces<aod::StoredJDielectronCollisionIds> storedDielectronCollisionIdsTable;
      Produces<aod::StoredDielectrons> storedDielectronsTable;
      Produces<aod::StoredJDielectronIds> storedDielectronIdsTable;
      Produces<aod::StoredDielectronsAll> storedDielectronsAllTable;
      Produces<aod::StoredJDielectronMcCollisions> storedDielectronMcCollisionsTable;
      Produces<aod::StoredJDielectronMcCollisionIds> storedDielectronMcCollisionIdsTable;
      Produces<aod::StoredJDielectronMcRCollDummys> storedDielectronMcRCollDummysTable;
      Produces<aod::StoredJDielectronMcs> storedDielectronParticlesTable;
      Produces<aod::StoredJDielectronMcIds> storedDielectronParticleIdsTable;
    } productsDielectron;

  } products;

  struct : PresliceGroup {
    Preslice<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>> TracksPerCollision = aod::jtrack::collisionId;

    Preslice<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
    Preslice<aod::McCollisionsD0> D0McCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsDplus> DplusMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsDs> DsMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsDstar> DstarMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsLc> LcMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsB0> B0McCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsBplus> BplusMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsXicToXiPiPi> XicToXiPiPiMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::McCollisionsDielectron> DielectronMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesD0MCP> D0ParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesDplusMCP> DplusParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesDsMCP> DsParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesDstarMCP> DstarParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesLcMCP> LcParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesB0MCP> B0ParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesBplusMCP> BplusParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    Preslice<aod::CandidatesXicToXiPiPiMCP> XicToXiPiPiParticlesPerMcCollision = aod::jcandidateindices::mcCollisionId;
    PresliceUnsorted<aod::JEMCTracks> EMCTrackPerTrack = aod::jemctrack::trackId;
  } preslices;

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
    if (config.performTrackSelection && !(track.trackSel() & ~((1 << jetderiveddatautilities::JTrackSel::trackSign) | (1 << jetderiveddatautilities::JTrackSel::notBadMcTrack)))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept. This includes also the track pT selction.
      return false;
    }
    if (track.pt() < config.trackPtSelectionMin || std::abs(track.eta()) > config.trackEtaSelectionMax) {
      return false;
    }
    return true;
  }

  template <bool isMc, typename T>
  void storeD0(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsD0 const& D0Collisions, T const& D0Candidates)
  {

    if (collision.isCollisionSelected()) {
      for (const auto& D0Collision : D0Collisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(D0Collision, products.productsD0.storedD0CollisionsTable);
        products.productsD0.storedD0CollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& D0Candidate : D0Candidates) {
        jethfutilities::fillHFCandidateTable<isMc>(D0Candidate, products.productsD0.storedD0CollisionsTable.lastIndex(), products.productsD0.storedD0sTable, products.productsD0.storedD0ParsTable, products.productsD0.storedD0ParExtrasTable, products.productsD0.storedD0ParDaughtersDummyTable, products.productsD0.storedD0SelsTable, products.productsD0.storedD0MlsTable, products.productsD0.storedD0MlDughtersDummyTable, products.productsD0.storedD0McsTable);
        products.productsD0.storedD0IdsTable(collisionMapping[collision.globalIndex()], trackMapping[D0Candidate.prong0Id()], trackMapping[D0Candidate.prong1Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeDplus(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDplus const& DplusCollisions, T const& DplusCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& DplusCollision : DplusCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(DplusCollision, products.productsDplus.storedDplusCollisionsTable);
        products.productsDplus.storedDplusCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& DplusCandidate : DplusCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(DplusCandidate, products.productsDplus.storedDplusCollisionsTable.lastIndex(), products.productsDplus.storedDplussTable, products.productsDplus.storedDplusParsTable, products.productsDplus.storedDplusParExtrasTable, products.productsDplus.storedDplusParDaughtersDummyTable, products.productsDplus.storedDplusSelsTable, products.productsDplus.storedDplusMlsTable, products.productsDplus.storedDplusMlDughtersDummyTable, products.productsDplus.storedDplusMcsTable);
        products.productsDplus.storedDplusIdsTable(collisionMapping[collision.globalIndex()], trackMapping[DplusCandidate.prong0Id()], trackMapping[DplusCandidate.prong1Id()], trackMapping[DplusCandidate.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeDs(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDs const& DsCollisions, T const& DsCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& DsCollision : DsCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(DsCollision, products.productsDs.storedDsCollisionsTable);
        products.productsDs.storedDsCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& DsCandidate : DsCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(DsCandidate, products.productsDs.storedDsCollisionsTable.lastIndex(), products.productsDs.storedDssTable, products.productsDs.storedDsParsTable, products.productsDs.storedDsParExtrasTable, products.productsDs.storedDsParDaughtersDummyTable, products.productsDs.storedDsSelsTable, products.productsDs.storedDsMlsTable, products.productsDs.storedDsMlDughtersDummyTable, products.productsDs.storedDsMcsTable);
        products.productsDs.storedDsIdsTable(collisionMapping[collision.globalIndex()], trackMapping[DsCandidate.prong0Id()], trackMapping[DsCandidate.prong1Id()], trackMapping[DsCandidate.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeDstar(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDstar const& DstarCollisions, T const& DstarCandidates)
  {

    if (collision.isCollisionSelected()) {
      for (const auto& DstarCollision : DstarCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(DstarCollision, products.productsDstar.storedDstarCollisionsTable);
        products.productsDstar.storedDstarCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& DstarCandidate : DstarCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(DstarCandidate, products.productsDstar.storedDstarCollisionsTable.lastIndex(), products.productsDstar.storedDstarsTable, products.productsDstar.storedDstarParsTable, products.productsDstar.storedDstarParExtrasDummyTable, products.productsDstar.storedDstarParD0sTable, products.productsDstar.storedDstarSelsTable, products.productsDstar.storedDstarMlsTable, products.productsDstar.storedDstarMlDaughtersDummyTable, products.productsDstar.storedDstarMcsTable);
        products.productsDstar.storedDstarIdsTable(collisionMapping[collision.globalIndex()], trackMapping[DstarCandidate.prong0Id()], trackMapping[DstarCandidate.prong1Id()], trackMapping[DstarCandidate.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeLc(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsLc const& LcCollisions, T const& LcCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& LcCollision : LcCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(LcCollision, products.productsLc.storedLcCollisionsTable);
        products.productsLc.storedLcCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& LcCandidate : LcCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(LcCandidate, products.productsLc.storedLcCollisionsTable.lastIndex(), products.productsLc.storedLcsTable, products.productsLc.storedLcParsTable, products.productsLc.storedLcParExtrasTable, products.productsLc.storedLcParDaughtersDummyTable, products.productsLc.storedLcSelsTable, products.productsLc.storedLcMlsTable, products.productsLc.storedLcMlDughtersDummyTable, products.productsLc.storedLcMcsTable);
        products.productsLc.storedLcIdsTable(collisionMapping[collision.globalIndex()], trackMapping[LcCandidate.prong0Id()], trackMapping[LcCandidate.prong1Id()], trackMapping[LcCandidate.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeB0(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsB0 const& B0Collisions, T const& B0Candidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& B0Collision : B0Collisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(B0Collision, products.productsB0.storedB0CollisionsTable);
        products.productsB0.storedB0CollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& B0Candidate : B0Candidates) {
        jethfutilities::fillHFCandidateTable<isMc>(B0Candidate, products.productsB0.storedB0CollisionsTable.lastIndex(), products.productsB0.storedB0sTable, products.productsB0.storedB0ParsTable, products.productsB0.storedB0ParExtrasTable, products.productsB0.storedB0ParDplussTable, products.productsB0.storedB0SelsTable, products.productsB0.storedB0MlsTable, products.productsB0.storedB0MlDplussTable, products.productsB0.storedB0McsTable);
        products.productsB0.storedB0IdsTable(collisionMapping[collision.globalIndex()], trackMapping[B0Candidate.prong0Id()], trackMapping[B0Candidate.prong1Id()], trackMapping[B0Candidate.prong2Id()], trackMapping[B0Candidate.prong3Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeBplus(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsBplus const& BplusCollisions, T const& BplusCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& BplusCollision : BplusCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(BplusCollision, products.productsBplus.storedBplusCollisionsTable);
        products.productsBplus.storedBplusCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& BplusCandidate : BplusCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(BplusCandidate, products.productsBplus.storedBplusCollisionsTable.lastIndex(), products.productsBplus.storedBplussTable, products.productsBplus.storedBplusParsTable, products.productsBplus.storedBplusParExtrasTable, products.productsBplus.storedBplusParD0sTable, products.productsBplus.storedBplusSelsTable, products.productsBplus.storedBplusMlsTable, products.productsBplus.storedBplusMlD0sTable, products.productsBplus.storedBplusMcsTable);
        products.productsBplus.storedBplusIdsTable(collisionMapping[collision.globalIndex()], trackMapping[BplusCandidate.prong0Id()], trackMapping[BplusCandidate.prong1Id()], trackMapping[BplusCandidate.prong2Id()]);
      }
    }
  }

  template <bool isMc, typename T>
  void storeXicToXiPiPi(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsXicToXiPiPi const& XicToXiPiPiCollisions, T const& XicToXiPiPiCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& XicToXiPiPiCollision : XicToXiPiPiCollisions) { // should only ever be one
        jethfutilities::fillHFCollisionTable(XicToXiPiPiCollision, products.productsXicToXiPiPi.storedXicToXiPiPiCollisionsTable);
        products.productsXicToXiPiPi.storedXicToXiPiPiCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& XicToXiPiPiCandidate : XicToXiPiPiCandidates) {
        jethfutilities::fillHFCandidateTable<isMc>(XicToXiPiPiCandidate, products.productsXicToXiPiPi.storedXicToXiPiPiCollisionsTable.lastIndex(), products.productsXicToXiPiPi.storedXicToXiPiPisTable, products.productsXicToXiPiPi.storedXicToXiPiPiParsTable, products.productsXicToXiPiPi.storedXicToXiPiPiParExtrasTable, products.productsXicToXiPiPi.storedXicToXiPiPiParDaughtersDummyTable, products.productsXicToXiPiPi.storedXicToXiPiPiSelsTable, products.productsXicToXiPiPi.storedXicToXiPiPiMlsTable, products.productsXicToXiPiPi.storedXicToXiPiPiMlDughtersDummyTable, products.productsXicToXiPiPi.storedXicToXiPiPiMcsTable);
        products.productsXicToXiPiPi.storedXicToXiPiPiIdsTable(collisionMapping[collision.globalIndex()], trackMapping[XicToXiPiPiCandidate.prong0Id()], trackMapping[XicToXiPiPiCandidate.prong1Id()], trackMapping[XicToXiPiPiCandidate.prong2Id()], trackMapping[XicToXiPiPiCandidate.prong3Id()], trackMapping[XicToXiPiPiCandidate.prong4Id()]);
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
  std::vector<int32_t> dsMcCollisionMapping;
  std::vector<int32_t> dstarMcCollisionMapping;
  std::vector<int32_t> lcMcCollisionMapping;
  std::vector<int32_t> b0McCollisionMapping;
  std::vector<int32_t> bplusMcCollisionMapping;
  std::vector<int32_t> xicToXiPiPiMcCollisionMapping;
  // std::vector<int32_t> dielectronMcCollisionMapping;

  void processBCs(soa::Join<aod::JCollisions, aod::JCollisionSelections> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const& bcs)
  {
    std::vector<int32_t> bcIndicies;
    bcMapping.clear();
    bcMapping.resize(bcs.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.triggerMask(), bc.timestamp(), bc.alias_raw(), bc.selection_raw(), bc.rct_raw());
          products.storedJBCParentIndexTable(bc.bcId());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping[bc.globalIndex()] = products.storedJBCsTable.lastIndex();
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBCs, "write out output tables for Bunch crossings", true);

  void processBCsForMcGenOnly(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::JBCs const& bcs)
  {
    std::vector<int32_t> bcIndicies;
    bcMapping.clear();
    bcMapping.resize(bcs.size(), -1);

    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        auto bc = mcCollision.bc_as<aod::JBCs>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.triggerMask(), bc.timestamp(), bc.alias_raw(), bc.selection_raw(), bc.rct_raw());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping[bc.globalIndex()] = products.storedJBCsTable.lastIndex();
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBCsForMcGenOnly, "write out output tables for Bunch crossings for mc gen only configurations", false);

  void processColllisons(soa::Join<aod::JCollisions, aod::JCollisionMcInfos, aod::JCollisionPIs, aod::JCollisionSelections> const& collisions)
  {
    collisionMapping.clear();
    collisionMapping.resize(collisions.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        products.storedJCollisionsTable(bcMapping[collision.bcId()], collision.posX(), collision.posY(), collision.posZ(), collision.collisionTime(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(), collision.centFV0A(), collision.centFV0M(), collision.centFT0A(), collision.centFT0C(), collision.centFT0M(), collision.centFT0CVariant1(), collision.hadronicRate(), collision.trackOccupancyInTimeRange(), collision.alias_raw(), collision.eventSel(), collision.rct_raw(), collision.triggerSel());
        collisionMapping[collision.globalIndex()] = products.storedJCollisionsTable.lastIndex();
        products.storedJCollisionMcInfosTable(collision.weight(), collision.getSubGeneratorId());
        products.storedJCollisionsParentIndexTable(collision.collisionId());
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processColllisons, "write out output tables for collisions", true);

  void processUPCCollisionInfos(soa::Join<aod::JCollisions, aod::JCollisionUPCs, aod::JCollisionSelections> const& collisions)
  {
    std::vector<float> amplitudesFV0;
    std::vector<float> amplitudesFT0A;
    std::vector<float> amplitudesFT0C;
    std::vector<float> amplitudesFDDA;
    std::vector<float> amplitudesFDDC;
    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        amplitudesFV0.clear();
        amplitudesFT0A.clear();
        amplitudesFT0C.clear();
        amplitudesFDDA.clear();
        amplitudesFDDC.clear();
        auto amplitudesFV0Span = collision.amplitudesFV0();
        auto amplitudesFT0ASpan = collision.amplitudesFT0A();
        auto amplitudesFT0CSpan = collision.amplitudesFT0C();
        auto amplitudesFDDASpan = collision.amplitudesFDDA();
        auto amplitudesFDDCSpan = collision.amplitudesFDDC();
        std::copy(amplitudesFV0Span.begin(), amplitudesFV0Span.end(), std::back_inserter(amplitudesFV0));
        std::copy(amplitudesFT0ASpan.begin(), amplitudesFT0ASpan.end(), std::back_inserter(amplitudesFT0A));
        std::copy(amplitudesFT0CSpan.begin(), amplitudesFT0CSpan.end(), std::back_inserter(amplitudesFT0C));
        std::copy(amplitudesFDDASpan.begin(), amplitudesFDDASpan.end(), std::back_inserter(amplitudesFDDA));
        std::copy(amplitudesFDDCSpan.begin(), amplitudesFDDCSpan.end(), std::back_inserter(amplitudesFDDC));
        products.storedJCollisionUPCsTable(amplitudesFV0, amplitudesFT0A, amplitudesFT0C, amplitudesFDDA, amplitudesFDDC);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processUPCCollisionInfos, "write out table for upc collision info", false);

  void processTracks(soa::Join<aod::JCollisions, aod::JCollisionSelections> const& collisions, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& tracks)
  {
    trackMapping.clear();
    trackMapping.resize(tracks.size(), -1);

    for (auto const& collision : collisions) {
      if (collision.isCollisionSelected()) {
        const auto tracksPerCollision = tracks.sliceBy(preslices.TracksPerCollision, collision.globalIndex());
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

  void processClusters(soa::Join<aod::JCollisions, aod::JEMCCollisionLbs, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::JEMCTracks const& emcTracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters)
  {
    if (collision.isCollisionSelected()) {
      products.storedJCollisionsEMCalLabelTable(collision.isAmbiguous(), collision.isEmcalReadout());
      for (const auto& cluster : clusters) {
        products.storedJClustersTable(collisionMapping[collision.globalIndex()], cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                      cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                      cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
        products.storedJClustersParentIndexTable(cluster.clusterId());

        std::vector<int32_t> clusterStoredJTrackIDs;
        for (const auto& clusterTrack : cluster.matchedTracks_as<aod::JTracks>()) {
          clusterStoredJTrackIDs.push_back(trackMapping[clusterTrack.globalIndex()]);
          auto emcTracksPerTrack = emcTracks.sliceBy(preslices.EMCTrackPerTrack, clusterTrack.globalIndex());
          auto emcTrackPerTrack = emcTracksPerTrack.iteratorAt(0);
          products.storedJTracksEMCalTable(trackMapping[clusterTrack.globalIndex()], emcTrackPerTrack.etaEmcal(), emcTrackPerTrack.phiEmcal(), emcTrackPerTrack.etaDiff(), emcTrackPerTrack.phiDiff());
        }
        products.storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processClusters, "write out output tables for clusters", false);

  //!!!!!!!!!! need to add the hadronic corrected energy and delete the new dummy process function

  void processD0Data(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0Data const& D0Candidates)
  {
    storeD0<false>(collision, tracks, D0Collisions, D0Candidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0Data, "write out data output tables for D0", false);

  void processD0MCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0MCD const& D0Candidates)
  {
    storeD0<true>(collision, tracks, D0Collisions, D0Candidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processD0MCD, "write out mcd output tables for D0", false);

  void processDplusData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDplus const& DplusCollisions, aod::CandidatesDplusData const& DplusCandidates)
  {
    storeDplus<false>(collision, tracks, DplusCollisions, DplusCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusData, "write out data output tables for Dplus", false);

  void processDplusMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDplus const& DplusCollisions, aod::CandidatesDplusMCD const& DplusCandidates)
  {
    storeDplus<true>(collision, tracks, DplusCollisions, DplusCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMCD, "write out mcd output tables for Dplus", false);

  void processDsData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDs const& DsCollisions, aod::CandidatesDsData const& DsCandidates)
  {
    storeDs<false>(collision, tracks, DsCollisions, DsCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDsData, "write out data output tables for Ds", false);

  void processDsMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDs const& DsCollisions, aod::CandidatesDsMCD const& DsCandidates)
  {
    storeDs<true>(collision, tracks, DsCollisions, DsCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDsMCD, "write out mcd output tables for Ds", false);

  void processDstarData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDstar const& DstarCollisions, aod::CandidatesDstarData const& DstarCandidates)
  {
    storeDstar<false>(collision, tracks, DstarCollisions, DstarCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDstarData, "write out data output tables for Dstar", false);

  void processDstarMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsDstar const& DstarCollisions, aod::CandidatesDstarMCD const& DstarCandidates)
  {
    storeDstar<true>(collision, tracks, DstarCollisions, DstarCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDstarMCD, "write out mcd output tables for Dstar", false);

  void processLcData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcData const& LcCandidates)
  {
    storeLc<false>(collision, tracks, LcCollisions, LcCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcData, "write out data output tables for Lc", false);

  void processLcMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcMCD const& LcCandidates)
  {
    storeLc<true>(collision, tracks, LcCollisions, LcCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcMCD, "write out mcd output tables for Lc", false);

  void processB0Data(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsB0 const& B0Collisions, aod::CandidatesB0Data const& B0Candidates)
  {
    storeB0<false>(collision, tracks, B0Collisions, B0Candidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processB0Data, "write out data output tables for B0", false);

  void processB0MCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsB0 const& B0Collisions, aod::CandidatesB0MCD const& B0Candidates)
  {
    storeB0<true>(collision, tracks, B0Collisions, B0Candidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processB0MCD, "write out mcd output tables for B0", false);

  void processBplusData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsBplus const& BplusCollisions, aod::CandidatesBplusData const& BplusCandidates)
  {
    storeBplus<false>(collision, tracks, BplusCollisions, BplusCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusData, "write out data output tables for bplus", false);

  void processBplusMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsBplus const& BplusCollisions, aod::CandidatesBplusMCD const& BplusCandidates)
  {
    storeBplus<true>(collision, tracks, BplusCollisions, BplusCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMCD, "write out mcd output tables for bplus", false);

  void processXicToXiPiPiData(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsXicToXiPiPi const& XicToXiPiPiCollisions, aod::CandidatesXicToXiPiPiData const& XicToXiPiPiCandidates)
  {
    storeXicToXiPiPi<false>(collision, tracks, XicToXiPiPiCollisions, XicToXiPiPiCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processXicToXiPiPiData, "write out data output tables for XicToXiPiPi", false);

  void processXicToXiPiPiMCD(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const& tracks, aod::CollisionsXicToXiPiPi const& XicToXiPiPiCollisions, aod::CandidatesXicToXiPiPiMCD const& XicToXiPiPiCandidates)
  {
    storeXicToXiPiPi<true>(collision, tracks, XicToXiPiPiCollisions, XicToXiPiPiCandidates);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processXicToXiPiPiMCD, "write out mcd output tables for XicToXiPiPi", false);

  void processDielectron(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, aod::JTracks const&, aod::CollisionsDielectron const& DielectronCollisions, aod::CandidatesDielectronData const& DielectronCandidates)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& DielectronCollision : DielectronCollisions) { // should only ever be one
        jetdqutilities::fillDielectronCollisionTable(DielectronCollision, products.productsDielectron.storedDielectronCollisionsTable);
        products.productsDielectron.storedDielectronCollisionIdsTable(collisionMapping[collision.globalIndex()]);
      }
      for (const auto& DielectronCandidate : DielectronCandidates) {
        jetdqutilities::fillDielectronCandidateTable(DielectronCandidate, products.productsDielectron.storedDielectronCollisionsTable.lastIndex(), products.productsDielectron.storedDielectronsTable, products.productsDielectron.storedDielectronsAllTable);
        products.productsDielectron.storedDielectronIdsTable(collisionMapping[collision.globalIndex()], trackMapping[DielectronCandidate.prong0Id()], trackMapping[DielectronCandidate.prong1Id()]);
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
        products.storedJMcCollisionsTable(bcMapping[mcCollision.bcId()], mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.multFV0A(), mcCollision.multFT0A(), mcCollision.multFT0C(), mcCollision.centFT0M(), mcCollision.weight(), mcCollision.accepted(), mcCollision.attempted(), mcCollision.xsectGen(), mcCollision.xsectErr(), mcCollision.ptHard(), mcCollision.rct_raw(), mcCollision.getGeneratorId(), mcCollision.getSubGeneratorId(), mcCollision.getSourceId(), mcCollision.impactParameter(), mcCollision.eventPlaneAngle());
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

        const auto particlesPerMcCollision = particles.sliceBy(preslices.ParticlesPerMcCollision, mcCollision.globalIndex());

        for (auto particle : particlesPerMcCollision) {
          particleMapping[particle.globalIndex()] = particleTableIndex;
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {

          std::vector<int32_t> mothersIds;
          int daughtersIds[2] = {-1, -1};
          if (config.savePartonLevelInfo) {
            if (particle.has_mothers()) {
              auto mothersIdTemps = particle.mothersIds();
              for (auto mothersIdTemp : mothersIdTemps) {
                mothersIds.push_back(particleMapping[mothersIdTemp]);
              }
            }
            if (particle.has_daughters()) {
              auto i = 0;
              for (auto daughterId : particle.daughtersIds()) {
                if (i > 1) {
                  break;
                }
                daughtersIds[i] = particleMapping[daughterId];
                i++;
              }
            }
          } else {
            if (!particle.isPhysicalPrimary()) { // add outgoing partons exclusion here later
              continue;
            }
          }
          products.storedJMcParticlesTable(mcCollisionMapping[mcCollision.globalIndex()], o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.statusCode(), particle.flags(), mothersIds, daughtersIds);
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
        const auto d0McCollisionsPerMcCollision = D0McCollisions.sliceBy(preslices.D0McCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& d0McCollisionPerMcCollision : d0McCollisionsPerMcCollision) {
          jethfutilities::fillHFMcCollisionTable(d0McCollisionPerMcCollision, products.productsD0.storedD0McCollisionsTable);
          products.productsD0.storedD0McCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          d0McCollisionMapping[d0McCollisionPerMcCollision.globalIndex()] = products.productsD0.storedD0McCollisionsTable.lastIndex();
        }
        const auto d0ParticlesPerMcCollision = D0Particles.sliceBy(preslices.D0ParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& D0Particle : d0ParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(D0Particle, products.productsD0.storedD0McCollisionsTable.lastIndex(), products.productsD0.storedD0ParticlesTable);
          products.productsD0.storedD0ParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[D0Particle.mcParticleId()]);
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
        const auto dplusMcCollisionsPerMcCollision = DplusMcCollisions.sliceBy(preslices.DplusMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& dplusMcCollisionPerMcCollision : dplusMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(dplusMcCollisionPerMcCollision, products.productsDplus.storedDplusMcCollisionsTable);
          products.productsDplus.storedDplusMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          dplusMcCollisionMapping[dplusMcCollisionPerMcCollision.globalIndex()] = products.productsDplus.storedDplusMcCollisionsTable.lastIndex();
        }
        const auto dplusParticlesPerMcCollision = DplusParticles.sliceBy(preslices.DplusParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& DplusParticle : dplusParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(DplusParticle, products.productsDplus.storedDplusMcCollisionsTable.lastIndex(), products.productsDplus.storedDplusParticlesTable);
          products.productsDplus.storedDplusParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DplusParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMCP, "write out Dplus mcp output tables", false);

  void processDsMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsDs const& DsMcCollisions, aod::CandidatesDsMCP const& DsParticles)
  {
    dsMcCollisionMapping.clear();
    dsMcCollisionMapping.resize(DsMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto dsMcCollisionsPerMcCollision = DsMcCollisions.sliceBy(preslices.DsMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& dsMcCollisionPerMcCollision : dsMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(dsMcCollisionPerMcCollision, products.productsDs.storedDsMcCollisionsTable);
          products.productsDs.storedDsMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          dsMcCollisionMapping[dsMcCollisionPerMcCollision.globalIndex()] = products.productsDs.storedDsMcCollisionsTable.lastIndex();
        }
        const auto dsParticlesPerMcCollision = DsParticles.sliceBy(preslices.DsParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& DsParticle : dsParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(DsParticle, products.productsDs.storedDsMcCollisionsTable.lastIndex(), products.productsDs.storedDsParticlesTable);
          products.productsDs.storedDsParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DsParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDsMCP, "write out Ds mcp output tables", false);

  void processDstarMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsDstar const& DstarMcCollisions, aod::CandidatesDstarMCP const& DstarParticles)
  {
    dstarMcCollisionMapping.clear();
    dstarMcCollisionMapping.resize(DstarMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto dstarMcCollisionsPerMcCollision = DstarMcCollisions.sliceBy(preslices.DstarMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& dstarMcCollisionPerMcCollision : dstarMcCollisionsPerMcCollision) {
          jethfutilities::fillHFMcCollisionTable(dstarMcCollisionPerMcCollision, products.productsDstar.storedDstarMcCollisionsTable);
          products.productsDstar.storedDstarMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          dstarMcCollisionMapping[dstarMcCollisionPerMcCollision.globalIndex()] = products.productsDstar.storedDstarMcCollisionsTable.lastIndex();
        }
        const auto dstarParticlesPerMcCollision = DstarParticles.sliceBy(preslices.DstarParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& DstarParticle : dstarParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(DstarParticle, products.productsDstar.storedDstarMcCollisionsTable.lastIndex(), products.productsDstar.storedDstarParticlesTable);
          products.productsDstar.storedDstarParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DstarParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDstarMCP, "write out D* mcp output tables", false);

  void processLcMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsLc const& LcMcCollisions, aod::CandidatesLcMCP const& LcParticles)
  {
    lcMcCollisionMapping.clear();
    lcMcCollisionMapping.resize(LcMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto lcMcCollisionsPerMcCollision = LcMcCollisions.sliceBy(preslices.LcMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& lcMcCollisionPerMcCollision : lcMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(lcMcCollisionPerMcCollision, products.productsLc.storedLcMcCollisionsTable);
          products.productsLc.storedLcMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          lcMcCollisionMapping[lcMcCollisionPerMcCollision.globalIndex()] = products.productsLc.storedLcMcCollisionsTable.lastIndex();
        }
        const auto lcParticlesPerMcCollision = LcParticles.sliceBy(preslices.LcParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& LcParticle : lcParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(LcParticle, products.productsLc.storedLcMcCollisionsTable.lastIndex(), products.productsLc.storedLcParticlesTable);
          products.productsLc.storedLcParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[LcParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcMCP, "write out Lc mcp output tables", false);

  void processB0MCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsB0 const& B0McCollisions, aod::CandidatesB0MCP const& B0Particles)
  {
    b0McCollisionMapping.clear();
    b0McCollisionMapping.resize(B0McCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto b0McCollisionsPerMcCollision = B0McCollisions.sliceBy(preslices.B0McCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& b0McCollisionPerMcCollision : b0McCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(b0McCollisionPerMcCollision, products.productsB0.storedB0McCollisionsTable);
          products.productsB0.storedB0McCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          b0McCollisionMapping[b0McCollisionPerMcCollision.globalIndex()] = products.productsB0.storedB0McCollisionsTable.lastIndex();
        }
        const auto b0ParticlesPerMcCollision = B0Particles.sliceBy(preslices.B0ParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& B0Particle : b0ParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(B0Particle, products.productsB0.storedB0McCollisionsTable.lastIndex(), products.productsB0.storedB0ParticlesTable);
          products.productsB0.storedB0ParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[B0Particle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processB0MCP, "write out B0 mcp output tables", false);

  void processBplusMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsBplus const& BplusMcCollisions, aod::CandidatesBplusMCP const& BplusParticles)
  {
    bplusMcCollisionMapping.clear();
    bplusMcCollisionMapping.resize(BplusMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto bplusMcCollisionsPerMcCollision = BplusMcCollisions.sliceBy(preslices.BplusMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& bplusMcCollisionPerMcCollision : bplusMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(bplusMcCollisionPerMcCollision, products.productsBplus.storedBplusMcCollisionsTable);
          products.productsBplus.storedBplusMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          bplusMcCollisionMapping[bplusMcCollisionPerMcCollision.globalIndex()] = products.productsBplus.storedBplusMcCollisionsTable.lastIndex();
        }
        const auto bplusParticlesPerMcCollision = BplusParticles.sliceBy(preslices.BplusParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& BplusParticle : bplusParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(BplusParticle, products.productsBplus.storedBplusMcCollisionsTable.lastIndex(), products.productsBplus.storedBplusParticlesTable);
          products.productsBplus.storedBplusParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[BplusParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMCP, "write out Bplus mcp output tables", false);

  void processXicToXiPiPiMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections> const& mcCollisions, aod::McCollisionsXicToXiPiPi const& XicToXiPiPiMcCollisions, aod::CandidatesXicToXiPiPiMCP const& XicToXiPiPiParticles)
  {
    xicToXiPiPiMcCollisionMapping.clear();
    xicToXiPiPiMcCollisionMapping.resize(XicToXiPiPiMcCollisions.size(), -1);
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isMcCollisionSelected()) {
        const auto xicToXiPiPiMcCollisionsPerMcCollision = XicToXiPiPiMcCollisions.sliceBy(preslices.XicToXiPiPiMcCollisionsPerMcCollision, mcCollision.globalIndex());
        for (const auto& xicToXiPiPiMcCollisionPerMcCollision : xicToXiPiPiMcCollisionsPerMcCollision) { // should only ever be one
          jethfutilities::fillHFMcCollisionTable(xicToXiPiPiMcCollisionPerMcCollision, products.productsXicToXiPiPi.storedXicToXiPiPiMcCollisionsTable);
          products.productsXicToXiPiPi.storedXicToXiPiPiMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
          xicToXiPiPiMcCollisionMapping[xicToXiPiPiMcCollisionPerMcCollision.globalIndex()] = products.productsXicToXiPiPi.storedXicToXiPiPiMcCollisionsTable.lastIndex();
        }
        const auto xicToXiPiPiParticlesPerMcCollision = XicToXiPiPiParticles.sliceBy(preslices.XicToXiPiPiParticlesPerMcCollision, mcCollision.globalIndex());
        for (const auto& XicToXiPiPiParticle : xicToXiPiPiParticlesPerMcCollision) {
          jethfutilities::fillHFCandidateMcTable(XicToXiPiPiParticle, products.productsXicToXiPiPi.storedXicToXiPiPiMcCollisionsTable.lastIndex(), products.productsXicToXiPiPi.storedXicToXiPiPiParticlesTable);
          products.productsXicToXiPiPi.storedXicToXiPiPiParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[XicToXiPiPiParticle.mcParticleId()]);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processXicToXiPiPiMCP, "write out XicToXiPiPi mcp output tables", false);

  void processDielectronMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, aod::JMcParticles const&, soa::Join<aod::McCollisionsDielectron, aod::JDielectronMcRCollDummys> const& DielectronMcCollisions, aod::CandidatesDielectronMCP const& DielectronParticles)
  {
    if (mcCollision.isMcCollisionSelected()) {

      const auto dielectronMcCollisionsPerMcCollision = DielectronMcCollisions.sliceBy(preslices.DielectronMcCollisionsPerMcCollision, mcCollision.globalIndex());
      for (const auto& dielectronMcCollisionPerMcCollision : dielectronMcCollisionsPerMcCollision) { // should only ever be one
        jetdqutilities::fillDielectronMcCollisionTable(dielectronMcCollisionPerMcCollision, products.productsDielectron.storedDielectronMcCollisionsTable);
        products.productsDielectron.storedDielectronMcCollisionIdsTable(mcCollisionMapping[mcCollision.globalIndex()]);
        products.productsDielectron.storedDielectronMcRCollDummysTable(dielectronMcCollisionPerMcCollision.dummyDQ());
      }
      for (const auto& DielectronParticle : DielectronParticles) {
        jetdqutilities::fillDielectronCandidateMcTable(DielectronParticle, products.productsDielectron.storedDielectronMcCollisionsTable.lastIndex(), products.productsDielectron.storedDielectronParticlesTable);
        std::vector<int32_t> DielectronMothersIds;
        int DielectronDaughtersIds[2] = {-1, -1};
        if (DielectronParticle.has_mothers()) {
          for (auto const& DielectronMother : DielectronParticle.template mothers_as<aod::JMcParticles>()) {
            DielectronMothersIds.push_back(particleMapping[DielectronMother.globalIndex()]);
          }
        }
        if (DielectronParticle.has_daughters()) {
          auto i = 0;
          for (auto const& DielectronDaughter : DielectronParticle.template daughters_as<aod::JMcParticles>()) {
            if (i > 1) {
              break;
            }
            DielectronDaughtersIds[i] = particleMapping[DielectronDaughter.globalIndex()];
            i++;
          }
        }
        products.productsDielectron.storedDielectronParticleIdsTable(mcCollisionMapping[mcCollision.globalIndex()], particleMapping[DielectronParticle.mcParticleId()], DielectronMothersIds, DielectronDaughtersIds);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDielectronMCP, "write out Dielectron mcp output tables", false);

  void processColllisonsMcCollisionLabel(soa::Join<aod::JCollisions, aod::JMcCollisionLbs, aod::JCollisionSelections>::iterator const& collision)
  {
    if (collision.isCollisionSelected()) {
      if (collision.has_mcCollision()) {
        products.storedJMcCollisionsLabelTable(mcCollisionMapping[collision.mcCollisionId()]);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processColllisonsMcCollisionLabel, "write out collision mcCollision label output tables", false);

  void processTracksMcParticleLabel(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, soa::Join<aod::JTracks, aod::JMcTrackLbs> const& tracks, aod::JMcParticles const&)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& track : tracks) {
        if (!trackSelection(track)) {
          continue;
        }
        if (track.has_mcParticle() && (config.savePartonLevelInfo || track.mcParticle().isPhysicalPrimary())) {
          products.storedJMcTracksLabelTable(particleMapping[track.mcParticleId()]);
        } else {
          products.storedJMcTracksLabelTable(-1);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processTracksMcParticleLabel, "write out track mcParticle label output tables", false);

  void processClusterMcLabel(soa::Join<aod::JCollisions, aod::JCollisionSelections>::iterator const& collision, soa::Join<aod::JClusters, aod::JMcClusterLbs> const& clusters, aod::JMcParticles const& particles)
  {
    if (collision.isCollisionSelected()) {
      for (const auto& cluster : clusters) {
        std::vector<int32_t> clusterStoredJParticleIDs;
        for (const auto& clusterParticleId : cluster.mcParticlesIds()) {
          if (!config.savePartonLevelInfo) {
            const auto& particle = particles.iteratorAt(clusterParticleId);
            if (!particle.isPhysicalPrimary()) {
              continue;
            }
          }
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
        products.productsD0.storedD0McCollisionsMatchingTable(d0CollisionIDs);
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
        products.productsDplus.storedDplusMcCollisionsMatchingTable(dplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDplusMcCollisionMatch, "write out Dplus McCollision collision label output tables", false);

  void processDsMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsDs, aod::HfDsMcRCollIds> const& DsMcCollisions, aod::CollisionsDs const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& DsMcCollision : DsMcCollisions) { // should just be one
        std::vector<int32_t> dplusCollisionIDs;
        for (auto const& dplusCollisionPerMcCollision : DsMcCollision.hfCollBases_as<aod::CollisionsDs>()) {
          dplusCollisionIDs.push_back(dplusMcCollisionMapping[dplusCollisionPerMcCollision.globalIndex()]);
        }
        products.productsDs.storedDsMcCollisionsMatchingTable(dplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDsMcCollisionMatch, "write out Ds McCollision collision label output tables", false);

  void processDstarMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsDstar, aod::HfDstarMcRCollIds> const& DstarMcCollisions, aod::CollisionsDstar const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& DstarMcCollision : DstarMcCollisions) { // should just be one
        std::vector<int32_t> dstarCollisionIDs;
        for (auto const& dstarCollisionPerMcCollision : DstarMcCollision.hfCollBases_as<aod::CollisionsDstar>()) {
          dstarCollisionIDs.push_back(dstarMcCollisionMapping[dstarCollisionPerMcCollision.globalIndex()]);
        }
        products.productsDstar.storedDstarMcCollisionsMatchingTable(dstarCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDstarMcCollisionMatch, "write out D* McCollision collision label output tables", false);

  void processLcMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsLc, aod::HfLcMcRCollIds> const& LcMcCollisions, aod::CollisionsLc const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& LcMcCollision : LcMcCollisions) { // should just be one
        std::vector<int32_t> lcCollisionIDs;
        for (auto const& lcCollisionPerMcCollision : LcMcCollision.hfCollBases_as<aod::CollisionsLc>()) {
          lcCollisionIDs.push_back(lcMcCollisionMapping[lcCollisionPerMcCollision.globalIndex()]);
        }
        products.productsLc.storedLcMcCollisionsMatchingTable(lcCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processLcMcCollisionMatch, "write out Lc McCollision collision label output tables", false);

  void processB0McCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsB0, aod::HfB0McRCollIds> const& B0McCollisions, aod::CollisionsB0 const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& B0McCollision : B0McCollisions) { // should just be one
        std::vector<int32_t> b0CollisionIDs;
        for (auto const& b0CollisionPerMcCollision : B0McCollision.hfCollBases_as<aod::CollisionsB0>()) {
          b0CollisionIDs.push_back(b0McCollisionMapping[b0CollisionPerMcCollision.globalIndex()]);
        }
        products.productsB0.storedB0McCollisionsMatchingTable(b0CollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processB0McCollisionMatch, "write out B0 McCollision collision label output tables", false);

  void processBplusMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsBplus, aod::HfBplusMcRCollIds> const& BplusMcCollisions, aod::CollisionsBplus const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& BplusMcCollision : BplusMcCollisions) { // should just be one
        std::vector<int32_t> bplusCollisionIDs;
        for (auto const& bplusCollisionPerMcCollision : BplusMcCollision.hfCollBases_as<aod::CollisionsBplus>()) {
          bplusCollisionIDs.push_back(bplusMcCollisionMapping[bplusCollisionPerMcCollision.globalIndex()]);
        }
        products.productsBplus.storedBplusMcCollisionsMatchingTable(bplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processBplusMcCollisionMatch, "write out Bplus McCollision collision label output tables", false);

  void processXicToXiPiPiMcCollisionMatch(soa::Join<aod::JMcCollisions, aod::JMcCollisionSelections>::iterator const& mcCollision, soa::Join<aod::McCollisionsXicToXiPiPi, aod::HfXicToXiPiPiMcRCollIds> const& XicToXiPiPiMcCollisions, aod::CollisionsXicToXiPiPi const&)
  {
    if (mcCollision.isMcCollisionSelected()) {
      for (const auto& XicToXiPiPiMcCollision : XicToXiPiPiMcCollisions) { // should just be one
        std::vector<int32_t> dplusCollisionIDs;
        for (auto const& dplusCollisionPerMcCollision : XicToXiPiPiMcCollision.hfCollBases_as<aod::CollisionsXicToXiPiPi>()) {
          dplusCollisionIDs.push_back(dplusMcCollisionMapping[dplusCollisionPerMcCollision.globalIndex()]);
        }
        products.productsXicToXiPiPi.storedXicToXiPiPiMcCollisionsMatchingTable(dplusCollisionIDs);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processXicToXiPiPiMcCollisionMatch, "write out XicToXiPiPi McCollision collision label output tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
