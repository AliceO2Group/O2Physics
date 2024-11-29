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

// task to produce a self contained data format for jet analyses from the full AO2D
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexer.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "EventFiltering/filterTables.h"
#include "EventFiltering/Zorro.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetV0Utilities.h"
#include "PWGJE/Core/JetDQUtilities.h"

#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataProducerTask {
  Produces<aod::BCCounts> bcCountsTable;
  Produces<aod::CollisionCounts> collisionCountsTable;
  Produces<aod::JDummys> jDummysTable;
  Produces<aod::JBCs> jBCsTable;
  Produces<aod::JBCPIs> jBCParentIndexTable;
  Produces<aod::JCollisions> jCollisionsTable;
  Produces<aod::JCollisionPIs> jCollisionsParentIndexTable;
  Produces<aod::JCollisionBCs> jCollisionsBunchCrossingIndexTable;
  Produces<aod::JEMCCollisionLbs> jCollisionsEMCalLabelTable;
  Produces<aod::JMcCollisionLbs> jMcCollisionsLabelTable;
  Produces<aod::JMcCollisions> jMcCollisionsTable;
  Produces<aod::JMcCollisionPIs> jMcCollisionsParentIndexTable;
  Produces<aod::JTracks> jTracksTable;
  Produces<aod::JTrackExtras> jTracksExtraTable;
  Produces<aod::JEMCTracks> jTracksEMCalTable;
  Produces<aod::JTrackPIs> jTracksParentIndexTable;
  Produces<aod::JMcTrackLbs> jMcTracksLabelTable;
  Produces<aod::JMcParticles> jMcParticlesTable;
  Produces<aod::JMcParticlePIs> jParticlesParentIndexTable;
  Produces<aod::JClusters> jClustersTable;
  Produces<aod::JClusterPIs> jClustersParentIndexTable;
  Produces<aod::JClusterTracks> jClustersMatchedTracksTable;
  Produces<aod::JMcClusterLbs> jMcClustersLabelTable;
  Produces<aod::JD0CollisionIds> jD0CollisionIdsTable;
  Produces<aod::JD0McCollisionIds> jD0McCollisionIdsTable;
  Produces<aod::JD0Ids> jD0IdsTable;
  Produces<aod::JD0PIds> jD0ParticleIdsTable;
  Produces<aod::JLcCollisionIds> jLcCollisionIdsTable;
  Produces<aod::JLcMcCollisionIds> jLcMcCollisionIdsTable;
  Produces<aod::JLcIds> jLcIdsTable;
  Produces<aod::JLcPIds> jLcParticleIdsTable;
  Produces<aod::JV0Ids> jV0IdsTable;
  Produces<aod::JV0McCollisions> jV0McCollisionsTable;
  Produces<aod::JV0McCollisionIds> jV0McCollisionIdsTable;
  Produces<aod::JV0Mcs> jV0McsTable;
  Produces<aod::JV0McIds> jV0McIdsTable;
  Produces<aod::JDielectronCollisionIds> jDielectronCollisionIdsTable;
  Produces<aod::JDielectronIds> jDielectronIdsTable;
  Produces<aod::JDielectronMcCollisions> jDielectronMcCollisionsTable;
  Produces<aod::JDielectronMcCollisionIds> jDielectronMcCollisionIdsTable;
  Produces<aod::JDielectronMcs> jDielectronMcsTable;
  Produces<aod::JDielectronMcIds> jDielectronMcIdsTable;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<float> dcaZMax{"dcaZMax", 0.2, "maximum DCAZ selection for tracks - only applied for reassociation"};

  Configurable<std::string> ccdbURL{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> includeTriggers{"includeTriggers", false, "fill the collision information with software trigger decisions"};

  Preslice<aod::EMCALClusterCells> perClusterCells = aod::emcalclustercell::emcalclusterId;
  Preslice<aod::EMCALMatchedTracks> perClusterTracks = aod::emcalclustercell::emcalclusterId;
  Preslice<aod::TrackAssoc> perCollisionTrackIndices = aod::track_association::collisionId;

  std::map<std::pair<int32_t, int32_t>, int32_t> trackCollisionMapping;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  Zorro triggerDecider;

  int runNumber;
  bool withCollisionAssociator;
  void init(InitContext const&)
  {
    if (doprocessTracksWithCollisionAssociator || includeTriggers) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      runNumber = 0;
      if (doprocessTracksWithCollisionAssociator) {
        withCollisionAssociator = true;
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
      } else {
        withCollisionAssociator = false;
      }
    }
  }

  void processClearMaps(aod::Collisions const&)
  {
    trackCollisionMapping.clear();
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processClearMaps, "clears all maps", true);

  void processBunchCrossings(soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>::iterator const& bc)
  {
    jBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp(), bc.alias_raw(), bc.selection_raw());
    jBCParentIndexTable(bc.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processBunchCrossings, "produces derived bunch crossing table", false);

  void processCollisions(soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::CentFT0Cs>::iterator const& collision, soa::Join<aod::BCs, aod::Timestamps> const&)
  {
    uint64_t triggerBit = 0;
    if (includeTriggers) {
      auto bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps>>();
      triggerDecider.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), jetderiveddatautilities::JTriggerMasks);
      triggerBit = jetderiveddatautilities::setTriggerSelectionBit(triggerDecider.getTriggerOfInterestResults(bc.globalBC()));
    }
    jCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0C(), collision.centFT0C(), collision.trackOccupancyInTimeRange(), jetderiveddatautilities::setEventSelectionBit(collision), collision.alias_raw(), triggerBit); // note change multFT0C to multFT0M when problems with multFT0A are fixed
    jCollisionsParentIndexTable(collision.globalIndex());
    jCollisionsBunchCrossingIndexTable(collision.bcId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processCollisions, "produces derived collision tables", true);

  void processCollisionsWithoutCentralityAndMultiplicity(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::BCs, aod::Timestamps> const&)
  {
    uint64_t triggerBit = 0;
    if (includeTriggers) {
      auto bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps>>();
      triggerDecider.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), jetderiveddatautilities::JTriggerMasks);
      triggerBit = jetderiveddatautilities::setTriggerSelectionBit(triggerDecider.getTriggerOfInterestResults(bc.globalBC()));
    }
    jCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), -1.0, -1.0, -1, jetderiveddatautilities::setEventSelectionBit(collision), collision.alias_raw(), triggerBit);
    jCollisionsParentIndexTable(collision.globalIndex());
    jCollisionsBunchCrossingIndexTable(collision.bcId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processCollisionsWithoutCentralityAndMultiplicity, "produces derived collision tables without centrality or multiplicity", false);

  void processCollisionsRun2(soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::CentRun2V0Ms>::iterator const& collision)
  {
    jCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0C(), collision.centRun2V0M(), -1, jetderiveddatautilities::setEventSelectionBit(collision), collision.alias_raw(), 0); // note change multFT0C to multFT0M when problems with multFT0A are fixed
    jCollisionsParentIndexTable(collision.globalIndex());
    jCollisionsBunchCrossingIndexTable(collision.bcId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processCollisionsRun2, "produces derived collision tables for Run 2 data", false);

  void processCollisionsALICE3(aod::Collision const& collision)
  {
    jCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), -1.0, -1.0, -1, -1.0, 0, 0);
    jCollisionsParentIndexTable(collision.globalIndex());
    jCollisionsBunchCrossingIndexTable(-1);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processCollisionsALICE3, "produces derived collision tables for ALICE 3 simulations", false);

  void processWithoutEMCalCollisionLabels(aod::Collision const&)
  {
    jCollisionsEMCalLabelTable(false, false);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processWithoutEMCalCollisionLabels, "produces dummy derived collision labels for EMCal", true);

  void processEMCalCollisionLabels(aod::EMCALMatchedCollision const& collision)
  {
    jCollisionsEMCalLabelTable(collision.ambiguous(), collision.isemcreadout());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processEMCalCollisionLabels, "produces derived collision labels for EMCal", false);

  void processMcCollisionLabels(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision)
  {

    if (collision.has_mcCollision()) {
      jMcCollisionsLabelTable(collision.mcCollisionId());
    } else {
      jMcCollisionsLabelTable(-1);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcCollisionLabels, "produces derived MC collision labels table", false);

  void processMcCollisions(aod::McCollision const& McCollision)
  {
    jMcCollisionsTable(McCollision.posX(), McCollision.posY(), McCollision.posZ(), McCollision.weight());
    jMcCollisionsParentIndexTable(McCollision.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcCollisions, "produces derived MC collision table", false);

  void processTracks(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::TrackSelectionExtension>::iterator const& track, aod::Collisions const&)
  {
    jTracksTable(track.collisionId(), track.pt(), track.eta(), track.phi(), jetderiveddatautilities::setTrackSelectionBit(track, track.dcaZ(), dcaZMax));
    auto trackParCov = getTrackParCov(track);
    auto xyzTrack = trackParCov.getXYZGlo();
    float sigmaDCAXYZ2;
    float dcaXYZ = getDcaXYZ(track, &sigmaDCAXYZ2);
    float dcaX = -99.0;
    float dcaY = -99.0;
    if (track.collisionId() >= 0) {
      auto const& collision = track.collision_as<aod::Collisions>();
      dcaX = xyzTrack.X() - collision.posX();
      dcaY = xyzTrack.Y() - collision.posY();
    }

    jTracksExtraTable(dcaX, dcaY, track.dcaZ(), track.dcaXY(), dcaXYZ, std::sqrt(track.sigmaDcaZ2()), std::sqrt(track.sigmaDcaXY2()), std::sqrt(sigmaDCAXYZ2), track.sigma1Pt()); // why is this getSigmaZY
    jTracksParentIndexTable(track.globalIndex());
    trackCollisionMapping[{track.globalIndex(), track.collisionId()}] = jTracksTable.lastIndex();
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processTracks, "produces derived track table", true);

  void processTracksWithCollisionAssociator(aod::Collisions const& collisions, soa::Join<aod::BCs, aod::Timestamps> const&, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::TrackSelectionExtension> const&, aod::TrackAssoc const& assocCollisions)
  {
    for (auto const& collision : collisions) {
      auto collisionTrackIndices = assocCollisions.sliceBy(perCollisionTrackIndices, collision.globalIndex());
      for (auto const& collisionTrackIndex : collisionTrackIndices) {
        auto track = collisionTrackIndex.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TracksDCACov, aod::TrackSelection, aod::TrackSelectionExtension>>();
        auto trackParCov = getTrackParCov(track);
        if (track.collisionId() == collision.globalIndex()) {
          jTracksTable(collision.globalIndex(), track.pt(), track.eta(), track.phi(), jetderiveddatautilities::setTrackSelectionBit(track, track.dcaZ(), dcaZMax));
          jTracksParentIndexTable(track.globalIndex());
          auto xyzTrack = trackParCov.getXYZGlo();
          float sigmaDCAXYZ2;
          float dcaXYZ = getDcaXYZ(track, &sigmaDCAXYZ2);
          jTracksExtraTable(xyzTrack.X() - collision.posX(), xyzTrack.Y() - collision.posY(), track.dcaZ(), track.dcaXY(), dcaXYZ, std::sqrt(track.sigmaDcaZ2()), std::sqrt(track.sigmaDcaXY2()), std::sqrt(sigmaDCAXYZ2), track.sigma1Pt()); // why is this getSigmaZY
        } else {
          auto bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps>>();
          initCCDB(bc, runNumber, ccdb, doprocessCollisionsRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, doprocessCollisionsRun2);
          o2::dataformats::DCA dcaCovInfo;
          dcaCovInfo.set(-999., -999., -999., -999., -999.);
          o2::dataformats::VertexBase collisionInfo;
          collisionInfo.setPos({collision.posX(), collision.posY(), collision.posZ()});
          collisionInfo.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(collisionInfo, trackParCov, 2.f, noMatCorr, &dcaCovInfo);
          jTracksTable(collision.globalIndex(), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), jetderiveddatautilities::setTrackSelectionBit(track, dcaCovInfo.getZ(), dcaZMax)); // only qualitytracksWDCA are a reliable selection
          jTracksParentIndexTable(track.globalIndex());
          auto xyzTrack = trackParCov.getXYZGlo();
          float dcaXY = dcaCovInfo.getY();
          float dcaZ = dcaCovInfo.getZ();
          float dcaXYZ = std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ);
          float covYY = dcaCovInfo.getSigmaY2();
          float covZZ = dcaCovInfo.getSigmaZ2();
          float covYZ = dcaCovInfo.getSigmaYZ();
          float sigmaDCAXYZ;
          if (dcaXYZ < o2::constants::math::Almost0) {
            sigmaDCAXYZ = o2::constants::math::VeryBig; // Protection against division by zero
          } else {
            sigmaDCAXYZ = covYY * (2.f * dcaXY / dcaXYZ) * (2.f * dcaXY / dcaXYZ) + covZZ * (2.f * dcaZ / dcaXYZ) * (2.f * dcaZ / dcaXYZ) + 2.f * covYZ * (2.f * dcaXY / dcaXYZ) * (2.f * dcaZ / dcaXYZ);
          }
          jTracksExtraTable(xyzTrack.X() - collision.posX(), xyzTrack.Y() - collision.posY(), dcaZ, dcaXY, dcaXYZ, std::sqrt(covZZ), std::sqrt(covYY), std::sqrt(sigmaDCAXYZ), std::sqrt(trackParCov.getSigma1Pt2()));
        }
        trackCollisionMapping[{track.globalIndex(), collision.globalIndex()}] = jTracksTable.lastIndex();
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processTracksWithCollisionAssociator, "produces derived track table taking into account track-to-collision associations", false);

  void processMcTrackLabels(soa::Join<aod::Tracks, aod::McTrackLabels>::iterator const& track)
  {
    if (track.has_mcParticle()) {
      jMcTracksLabelTable(track.mcParticleId());
    } else {
      jMcTracksLabelTable(-1);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcTrackLabels, "produces derived track labels table", false);

  void processMcTrackLabelsWithCollisionAssociator(aod::Collisions const& collisions, soa::Join<aod::Tracks, aod::McTrackLabels> const&, aod::TrackAssoc const& assocCollisions)
  {
    for (auto const& collision : collisions) {
      auto collisionTrackIndices = assocCollisions.sliceBy(perCollisionTrackIndices, collision.globalIndex());
      for (auto const& collisionTrackIndex : collisionTrackIndices) {
        auto track = collisionTrackIndex.track_as<soa::Join<aod::Tracks, aod::McTrackLabels>>();
        if (track.collisionId() == collision.globalIndex() && track.has_mcParticle()) {
          jMcTracksLabelTable(track.mcParticleId());
        } else {
          jMcTracksLabelTable(-1);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcTrackLabelsWithCollisionAssociator, "produces derived track labels table taking into account track-to-collision associations", false);

  void processParticles(aod::McParticle const& particle)
  {
    std::vector<int32_t> mothersId;
    if (particle.has_mothers()) {
      auto mothersIdTemps = particle.mothersIds();
      for (auto mothersIdTemp : mothersIdTemps) {
        mothersId.push_back(mothersIdTemp);
      }
    }
    int daughtersId[2] = {-1, -1};
    auto i = 0;
    if (particle.has_daughters()) {
      for (auto daughterId : particle.daughtersIds()) {
        if (i > 1) {
          break;
        }
        daughtersId[i] = daughterId;
        i++;
      }
    }
    jMcParticlesTable(particle.mcCollisionId(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
    jParticlesParentIndexTable(particle.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processParticles, "produces derived parrticle table", false);

  void processClusters(aod::Collision const&, aod::EMCALClusters const& clusters, aod::EMCALClusterCells const& cells, aod::Calos const&, aod::EMCALMatchedTracks const& matchedTracks, soa::Join<aod::Tracks, aod::TracksExtra> const&)
  {

    for (auto cluster : clusters) {

      auto const clusterCells = cells.sliceBy(perClusterCells, cluster.globalIndex());

      float leadingCellEnergy = -1.0;
      float subleadingCellEnergy = -1.0;
      float cellAmplitude = -1.0;
      int leadingCellNumber = -1;
      int subleadingCellNumber = -1;
      int cellNumber = -1;
      for (auto const& clutserCell : clusterCells) {
        cellAmplitude = clutserCell.calo().amplitude();
        cellNumber = clutserCell.calo().cellNumber();
        if (cellAmplitude > subleadingCellEnergy) {
          subleadingCellEnergy = cellAmplitude;
          subleadingCellNumber = cellNumber;
        }
        if (subleadingCellEnergy > leadingCellEnergy) {
          std::swap(leadingCellEnergy, subleadingCellEnergy);
          std::swap(leadingCellNumber, subleadingCellNumber);
        }
      }

      jClustersTable(cluster.collisionId(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(), cluster.nlm(), cluster.definition(), leadingCellEnergy, subleadingCellEnergy, leadingCellNumber, subleadingCellNumber);
      jClustersParentIndexTable(cluster.globalIndex());

      auto const clusterTracks = matchedTracks.sliceBy(perClusterTracks, cluster.globalIndex());
      std::vector<int32_t> clusterTrackIDs;
      for (const auto& clusterTrack : clusterTracks) {
        auto JClusterID = trackCollisionMapping.find({clusterTrack.trackId(), cluster.collisionId()}); // does EMCal use its own associator?
        clusterTrackIDs.push_back(JClusterID->second);
        auto emcTrack = clusterTrack.track_as<soa::Join<aod::Tracks, aod::TracksExtra>>();
        jTracksEMCalTable(JClusterID->second, emcTrack.trackEtaEmcal(), emcTrack.trackPhiEmcal());
      }
      jClustersMatchedTracksTable(clusterTrackIDs);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processClusters, "produces derived cluster tables", false);

  void processMcClusterLabels(aod::EMCALMCCluster const& cluster)
  {
    std::vector<int32_t> particleIds;
    for (auto particleId : cluster.mcParticleIds()) {
      particleIds.push_back(particleId);
    }
    std::vector<float> amplitudeA;
    auto amplitudeASpan = cluster.amplitudeA();
    std::copy(amplitudeASpan.begin(), amplitudeASpan.end(), std::back_inserter(amplitudeA));
    jMcClustersLabelTable(particleIds, amplitudeA);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcClusterLabels, "produces derived cluster particle label table", false);

  void processD0Collisions(aod::HfD0CollIds::iterator const& D0Collision)
  {
    jD0CollisionIdsTable(D0Collision.collisionId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processD0Collisions, "produces derived index for D0 collisions", false);

  void processD0McCollisions(aod::HfD0McCollIds::iterator const& D0McCollision)
  {
    jD0McCollisionIdsTable(D0McCollision.mcCollisionId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processD0McCollisions, "produces derived index for D0 MC collisions", false);

  void processD0(aod::HfD0Ids::iterator const& D0, aod::Tracks const&)
  {
    auto JProng0ID = trackCollisionMapping.find({D0.prong0Id(), D0.prong0_as<aod::Tracks>().collisionId()});
    auto JProng1ID = trackCollisionMapping.find({D0.prong1Id(), D0.prong1_as<aod::Tracks>().collisionId()});
    if (withCollisionAssociator) {
      JProng0ID = trackCollisionMapping.find({D0.prong0Id(), D0.collisionId()});
      JProng1ID = trackCollisionMapping.find({D0.prong1Id(), D0.collisionId()});
    }
    jD0IdsTable(D0.collisionId(), JProng0ID->second, JProng1ID->second);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processD0, "produces derived index for D0 candidates", false);

  void processD0MC(aod::HfD0PIds::iterator const& D0)
  {
    jD0ParticleIdsTable(D0.mcCollisionId(), D0.mcParticleId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processD0MC, "produces derived index for D0 particles", false);

  void processLcCollisions(aod::Hf3PCollIds::iterator const& LcCollision)
  {
    jLcCollisionIdsTable(LcCollision.collisionId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processLcCollisions, "produces derived index for Lc collisions", false);

  void processLcMcCollisions(aod::Hf3PMcCollIds::iterator const& LcMcCollision)
  {
    jLcMcCollisionIdsTable(LcMcCollision.mcCollisionId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processLcMcCollisions, "produces derived index for Lc MC collisions", false);

  void processLc(aod::Hf3PIds::iterator const& Lc, aod::Tracks const&)
  {
    auto JProng0ID = trackCollisionMapping.find({Lc.prong0Id(), Lc.prong0_as<aod::Tracks>().collisionId()});
    auto JProng1ID = trackCollisionMapping.find({Lc.prong1Id(), Lc.prong1_as<aod::Tracks>().collisionId()});
    auto JProng2ID = trackCollisionMapping.find({Lc.prong2Id(), Lc.prong2_as<aod::Tracks>().collisionId()});
    if (withCollisionAssociator) {
      JProng0ID = trackCollisionMapping.find({Lc.prong0Id(), Lc.collisionId()});
      JProng1ID = trackCollisionMapping.find({Lc.prong1Id(), Lc.collisionId()});
      JProng2ID = trackCollisionMapping.find({Lc.prong2Id(), Lc.collisionId()});
    }
    jLcIdsTable(Lc.collisionId(), JProng0ID->second, JProng1ID->second, JProng2ID->second);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processLc, "produces derived index for Lc candidates", false);

  void processLcMC(aod::Hf3PPIds::iterator const& Lc)
  {
    jLcParticleIdsTable(Lc.mcCollisionId(), Lc.mcParticleId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processLcMC, "produces derived index for Lc particles", false);

  void processV0(aod::V0Indices::iterator const& V0, aod::Tracks const&)
  {
    auto JPosTrackID = trackCollisionMapping.find({V0.posTrackId(), V0.posTrack_as<aod::Tracks>().collisionId()});
    auto JNegTrackID = trackCollisionMapping.find({V0.negTrackId(), V0.negTrack_as<aod::Tracks>().collisionId()});
    if (withCollisionAssociator) {
      JPosTrackID = trackCollisionMapping.find({V0.posTrackId(), V0.collisionId()});
      JNegTrackID = trackCollisionMapping.find({V0.negTrackId(), V0.collisionId()});
    }
    jV0IdsTable(V0.collisionId(), JPosTrackID->second, JNegTrackID->second);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processV0, "produces derived index for V0 candidates", false);

  void processV0MC(aod::McCollision const& mcCollision, aod::McParticles const& particles)
  { // can loop over McV0Labels tables if we want to only store matched V0Particles
    bool filledV0McCollisionTable = false;
    for (auto const& particle : particles) {
      if (jetv0utilities::isV0Particle(particles, particle)) {
        if (!filledV0McCollisionTable) {
          jV0McCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
          jV0McCollisionIdsTable(mcCollision.globalIndex());
          filledV0McCollisionTable = true;
        }
        std::vector<int32_t> mothersId;
        if (particle.has_mothers()) {
          auto mothersIdTemps = particle.mothersIds();
          for (auto mothersIdTemp : mothersIdTemps) {
            mothersId.push_back(mothersIdTemp);
          }
        }
        int daughtersId[2] = {-1, -1};
        auto i = 0;
        if (particle.has_daughters()) {
          for (auto daughterId : particle.daughtersIds()) {
            if (i > 1) {
              break;
            }
            daughtersId[i] = daughterId;
            i++;
          }
        }
        auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
        jV0McsTable(jV0McCollisionsTable.lastIndex(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), pdgParticle->Mass(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), jetv0utilities::setV0ParticleDecayBit(particles, particle));
        jV0McIdsTable(mcCollision.globalIndex(), particle.globalIndex(), mothersId, daughtersId);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processV0MC, "produces V0 particles", false);

  void processDielectronCollisions(aod::ReducedEventsInfo::iterator const& DielectronCollision)
  {
    jDielectronCollisionIdsTable(DielectronCollision.collisionId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processDielectronCollisions, "produces derived index for Dielectron collisions", false);

  void processDielectron(aod::DielectronInfo const& Dielectron, aod::Tracks const&)
  {
    auto JProng0ID = trackCollisionMapping.find({Dielectron.prong0Id(), Dielectron.prong0_as<aod::Tracks>().collisionId()});
    auto JProng1ID = trackCollisionMapping.find({Dielectron.prong1Id(), Dielectron.prong1_as<aod::Tracks>().collisionId()});
    if (withCollisionAssociator) {
      JProng0ID = trackCollisionMapping.find({Dielectron.prong0Id(), Dielectron.collisionId()});
      JProng1ID = trackCollisionMapping.find({Dielectron.prong1Id(), Dielectron.collisionId()});
    }
    jDielectronIdsTable(Dielectron.collisionId(), JProng0ID->second, JProng1ID->second);
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processDielectron, "produces derived index for Dielectron candidates", false);

  void processDielectronMc(aod::McCollision const& mcCollision, aod::McParticles const& particles)
  {
    bool filledDielectronMcCollisionTable = false;
    for (auto const& particle : particles) {
      if (jetdqutilities::isDielectronParticle(particles, particle)) {
        if (!filledDielectronMcCollisionTable) {
          jDielectronMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
          jDielectronMcCollisionIdsTable(mcCollision.globalIndex());
          filledDielectronMcCollisionTable = true;
        }
        std::vector<int32_t> mothersId;
        if (particle.has_mothers()) {
          auto mothersIdTemps = particle.mothersIds();
          for (auto mothersIdTemp : mothersIdTemps) {
            mothersId.push_back(mothersIdTemp);
          }
        }
        int daughtersId[2] = {-1, -1};
        auto i = 0;
        if (particle.has_daughters()) {
          for (auto daughterId : particle.daughtersIds()) {
            if (i > 1) {
              break;
            }
            daughtersId[i] = daughterId;
            i++;
          }
        }
        auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
        jDielectronMcsTable(jDielectronMcCollisionsTable.lastIndex(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), pdgParticle->Mass(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), jetdqutilities::setDielectronParticleDecayBit(particles, particle), RecoDecay::getCharmHadronOrigin(particles, particle, false)); // Todo: should the last thing be false?
        jDielectronMcIdsTable(mcCollision.globalIndex(), particle.globalIndex(), mothersId, daughtersId);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processDielectronMc, "produces Dielectron mccollisions and particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerTask>(cfgc, TaskName{"jet-deriveddata-producer"})};
}
