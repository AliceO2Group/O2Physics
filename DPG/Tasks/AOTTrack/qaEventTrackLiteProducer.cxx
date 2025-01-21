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

///
/// \file   qaEventTrackLiteProducer.cxx
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Mattia Faggin <mattia.faggin@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief  Task to produce a table with reduced information used for correlation studies for track selection, ideally used with qaEventTrackite
///

#include "qaEventTrack.h"

#include <vector>

#include "TRandom.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct QaEventTrackLiteProducer {
  // Tables to produce
  Produces<o2::aod::DPGCollisions> tableCollisions;
  Produces<o2::aod::DPGCollsBig> tableCollsBig;
  Produces<o2::aod::DPGTracks> tableTracks;
  Produces<o2::aod::DPGRecoParticles> tableRecoParticles;
  Produces<o2::aod::DPGNonRecoParticles> tableNonRecoParticles;

  // general steering settings
  Configurable<bool> isRun3{"isRun3", false, "Is Run3 dataset"}; // TODO: derive this from metadata once possible to get rid of the flag

  // options to select specific events
  Configurable<bool> selectGoodEvents{"selectGoodEvents", true, "select good events"};
  // selection specific to the table creation workflow
  Configurable<float> selectMaxVtxZ{"selectMaxVtxZ", 100.f, "Derived data option: select collision in a given Z window"};
  Configurable<int> targetNumberOfEvents{"targetNumberOfEvents", 10000000, "Derived data option: target number of collisions, if the target is met, future collisions will be skipped"};
  Configurable<float> fractionOfSampledEvents{"fractionOfSampledEvents", 1.f, "Derived data option: fraction of events to sample"};
  Configurable<bool> storeOnlySinglePvCollsBig{"storeOnlySinglePvCollsBig", false, "Colls. big table (MC): store only reco. PV for single-reconstructed MC collisions"};
  Configurable<bool> storeOnlyMultiplePvCollsBig{"storeOnlyMultiplePvCollsBig", false, "Colls. big table (MC): store only reco. PV for multiple-reconstructed MC collisions"};

  // options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<int> selectCharge{"selectCharge", 0, "select charge +1 or -1 (0 means no selection)"};
  Configurable<bool> selectPrim{"selectPrim", false, "select primaries"};
  Configurable<bool> selectSec{"selectSec", false, "select secondaries"};
  Configurable<int> selectPID{"selectPID", 0, "select pid"};
  Configurable<float> minPt{"minPt", -10.f, "Minimum pt of accepted tracks"};
  Configurable<float> maxPt{"maxPt", 1e10f, "Maximum pt of accepted tracks"};
  Configurable<float> minEta{"minEta", -2.f, "Minimum eta of accepted tracks"};
  Configurable<float> maxEta{"maxEta", 2.0f, "Maximum eta of accepted tracks"};
  Configurable<float> minPhi{"minPhi", -1.f, "Minimum phi of accepted tracks"};
  Configurable<float> maxPhi{"maxPhi", 10.f, "Maximum phi of accepted tracks"};

  // TODO: ask if one can have different filters for both process functions
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  using TrackIUTable = soa::Join<aod::TracksIU, aod::TrackSelection>;
  Partition<TrackIUTable> tracksIUFiltered = (trackSelection.node() == 0) ||
                                             ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                                             ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                                             ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                                             ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                                             ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;
  Partition<TrackTableData> tracksFilteredCorrIU = (trackSelection.node() == 0) ||
                                                   ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                                                   ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                                                   ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                                                   ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                                                   ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<TrackTableData> perRecoCollision = aod::track::collisionId;

  int counterColl;
  int counterDF;

  void init(InitContext const&)
  {
    int howManyProcesses = static_cast<int>(doprocessTableData) + static_cast<int>(doprocessTableMC) + static_cast<int>(doprocessTableDataCollsBig) + static_cast<int>(doprocessTableMCCollsBig);
    if (howManyProcesses > 1) {
      LOGF(fatal, "%d process functions enabled. Enable only one of them!", howManyProcesses);
    } else if (howManyProcesses == 0) {
      LOGF(fatal, "No process function enabled. Enable one of them.");
    }

    /// for studies with collision table
    counterColl = 0;
    counterDF = 0;
    if (doprocessTableMCCollsBig && storeOnlySinglePvCollsBig && storeOnlyMultiplePvCollsBig) {
      LOGF(fatal, "storeOnlySinglePvCollsBig and storeOnlyMultiplePvCollsBig are both activated. Not possible. Fix the configuration.");
    }
  }

  // Function to select tracks
  template <bool IS_MC, typename T>
  bool isSelectedTrack(const T& track)
  {
    if (track.pt() < minPt || track.pt() > maxPt) { // Extra pT selection
      return false;
    }
    if (track.eta() < minEta || track.eta() > maxEta) { // Extra Eta selection
      return false;
    }
    if (track.phi() < minPhi || track.phi() > maxPhi) { // Extra Phi selection
      return false;
    }
    if (selectCharge && (selectCharge != track.sign())) {
      return false;
    }
    if constexpr (IS_MC) {
      if (!track.has_mcParticle()) {
        if (selectPrim || selectSec || selectPID) {
          return false;
        } else {
          return true;
        }
      }
      auto particle = track.mcParticle();
      const bool isPrimary = particle.isPhysicalPrimary();
      if (selectPrim && !isPrimary) {
        return false;
      }
      if (selectSec && isPrimary) {
        return false;
      }
      if (selectPID && selectPID != std::abs(particle.pdgCode())) {
        return false;
      }
    }
    return true;
  }

  // Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    if (selectGoodEvents && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    return true;
  }

  // Process function for data
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  using TrackTableMC = soa::Join<TrackTableData, aod::McTrackLabels>;

  // Process functions for skimming data
  void processTableData(CollisionTableData::iterator const& collision,
                        soa::Filtered<soa::Join<TrackTableData, aod::TOFSignal, aod::TOFEvTime>> const& tracks,
                        aod::BCs const& bcs)
  {
    fillDerivedTable<false>(collision, tracks, 0, bcs);
  }
  PROCESS_SWITCH(QaEventTrackLiteProducer, processTableData, "Process data for table producing", true);

  void processTableMC(CollisionTableMC::iterator const& collision,
                      soa::Filtered<soa::Join<TrackTableMC, aod::TOFSignal, aod::TOFEvTime>> const& tracks,
                      aod::McParticles const& mcParticles,
                      aod::McCollisions const&,
                      aod::BCs const& bcs)
  {
    fillDerivedTable<true>(collision, tracks, mcParticles, bcs);
  }
  PROCESS_SWITCH(QaEventTrackLiteProducer, processTableMC, "Process MC for table producing", false);

  //**************************************************************************************************
  /**
   * Fill reco level tables.
   */
  //**************************************************************************************************
  int nTableEventCounter = 0; // Number of processed events
  template <bool IS_MC, typename C, typename T, typename P>
  void fillDerivedTable(const C& collision, const T& tracks, const P& particles, const aod::BCs&)
  {
    if (!isSelectedCollision(collision)) {
      return;
    }
    if (std::abs(collision.posZ()) > selectMaxVtxZ) {
      return;
    }
    if (fractionOfSampledEvents < 1.f && (gRandom->Uniform()) > fractionOfSampledEvents) { // Skip events that are not sampled
      return;
    }
    if (nTableEventCounter > targetNumberOfEvents) { // Skip events if target is reached
      return;
    }
    nTableEventCounter++;

    tableCollisions(collision.posZ(),
                    (isRun3 ? collision.sel8() : collision.sel7()),
                    collision.bc().runNumber(), collision.numContrib());
    int nTracks = 0;
    int particleProduction = 0;

    for (const auto& track : tracks) {
      if (!isSelectedTrack<IS_MC>(track)) {
        continue;
      }
      ++nTracks;
    }
    tableTracks.reserve(nTracks);
    std::vector<int64_t> recoPartIndices(IS_MC ? nTracks : 0);

    if constexpr (IS_MC) { // Running only on MC
      tableRecoParticles.reserve(nTracks);
    }
    int64_t iTrack = 0;
    for (const auto& track : tracks) {
      if (!isSelectedTrack<IS_MC>(track)) {
        continue;
      }
      tableTracks(tableCollisions.lastIndex(),
                  track.pt(), track.tpcInnerParam(), track.eta(), track.phi(), track.pt() * std::sqrt(track.c1Pt21Pt2()),
                  track.flags(), track.sign(),
                  track.dcaXY(), track.dcaZ(), track.length(),
                  track.itsClusterMap(),
                  track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                  track.hasITS(), track.hasTPC(), track.hasTRD(), track.hasTOF(),
                  track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                  track.tpcCrossedRowsOverFindableCls(), track.tpcFoundOverFindableCls(), track.tpcFractionSharedCls(),
                  track.itsNCls(), track.itsNClsInnerBarrel(), track.tpcSignal(), track.tofSignal() - track.tofEvTime());

      if constexpr (IS_MC) { // Running only on MC
        if (track.has_mcParticle()) {
          auto particle = track.mcParticle();
          recoPartIndices[iTrack++] = particle.globalIndex();
          if (particle.isPhysicalPrimary()) {
            particleProduction = 0;
          } else if (particle.getProcess() == 4) {
            particleProduction = 1;
          } else {
            particleProduction = 2;
          }
          tableRecoParticles(particle.pt(), particle.eta(), particle.phi(),
                             particle.pdgCode(), particleProduction);
        } else { // If it does not have the particle we fill with the track values and tag it with -1 in the production
          tableRecoParticles(track.pt(), track.eta(), track.phi(),
                             0, -1);
        }
      }
    }

    // Running only on MC
    if constexpr (IS_MC) {
      if (!collision.has_mcCollision()) {
        return;
      }
      const auto& particlesInCollision = particles.sliceBy(perMcCollision, collision.mcCollision().globalIndex());
      tableNonRecoParticles.reserve(particlesInCollision.size() - nTracks);
      for (const auto& particle : particlesInCollision) {
        const auto partReconstructed = std::find(recoPartIndices.begin(), recoPartIndices.end(), particle.globalIndex()) != recoPartIndices.end();
        if (partReconstructed) {
          continue;
        }
        if (particle.isPhysicalPrimary()) {
          particleProduction = 0;
        } else if (particle.getProcess() == 4) {
          particleProduction = 1;
        } else {
          particleProduction = 2;
        }
        tableNonRecoParticles(tableCollisions.lastIndex(),
                              particle.pt(), particle.eta(), particle.phi(),
                              particle.pdgCode(), particleProduction,
                              particle.vx(), particle.vy(), particle.vz());
      }
    }
  }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  ///   Table filling for offline collision monitoring   ///
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
  using CollsBigTable = soa::Join<CollisionTableData, aod::Mults>;
  using CollsBigTableMC = soa::Join<CollsBigTable, aod::McCollisionLabels>;
  using McCollsWithExtra = soa::Join<aod::McCollisions, aod::McCollsExtra>;

  template <bool IS_MC, typename COLLS, typename MCCOLLS, typename TFILT, typename TALL>
  void fillCollsBigTable(COLLS& collisions, MCCOLLS&, TFILT& tracksFiltered, TALL& tracksAll, soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse> const&, aod::FT0s const&)
  {
    if (nTableEventCounter > targetNumberOfEvents) { // Skip events if target is reached
      return;
    }
    for (const auto& collision : collisions) {

      if (fractionOfSampledEvents < 1.f && (gRandom->Uniform()) > fractionOfSampledEvents) { // Skip events that are not sampled
        return;
      }
      nTableEventCounter++;

      const auto tracksFilteredPerColl = tracksFiltered.sliceBy(perRecoCollision, collision.globalIndex());
      const auto tracksAllPerColl = tracksAll.sliceBy(perRecoCollision, collision.globalIndex());

      const auto& bc = collision.template bc_as<BCsWithRun3Matchings>();
      int collIDMC = -1;
      float posXMC = -9999.;
      float posYMC = -9999.;
      float posZMC = -9999.;
      float collTimeMC = -9999.;
      int isFakeCollision = -1;
      int recoPVsPerMcColl = -1;            // how many reco. PV are associated to this MC collision
      int isPvHighestContribForMcColl = -1; // the "best" PV associated to this MC collision, i.e. the one with the highest number of PV contributors
      if constexpr (IS_MC) {
        if (!collision.has_mcCollision()) {
          isFakeCollision = 1;
        } else {
          isFakeCollision = 0;
          const auto& mcCollision = collision.template mcCollision_as<MCCOLLS>();
          collIDMC = mcCollision.globalIndex();
          posXMC = mcCollision.posX();
          posYMC = mcCollision.posY();
          posZMC = mcCollision.posZ();
          collTimeMC = mcCollision.t();
          recoPVsPerMcColl = mcCollision.numRecoCollision();
          isPvHighestContribForMcColl = (mcCollision.bestCollisionIndex() == collision.globalIndex());
        }

        /// MC
        /// Decide to fill the table only with the desired collisions
        if ((!storeOnlySinglePvCollsBig && !storeOnlyMultiplePvCollsBig) || /// store all collisions
            (storeOnlySinglePvCollsBig && recoPVsPerMcColl == 1) ||         /// store only single-reconstructed collisions
            (storeOnlyMultiplePvCollsBig && recoPVsPerMcColl > 1))          /// store only multiple-reconstructed collisions
        {
          tableCollsBig((isRun3 ? collision.sel8() : collision.sel7()),
                        bc.runNumber(),
                        collision.posX(), collision.posY(), collision.posZ(),
                        collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(),
                        collision.numContrib(), tracksAllPerColl.size(), tracksFilteredPerColl.size(),
                        collision.chi2(),
                        bc.globalBC(),
                        bc.has_ft0() ? bc.ft0().posZ() : -999.,
                        collision.multFT0A(), collision.multFT0C(), collision.multFT0M(), collision.multFV0A(),
                        collision.collisionTime(), collision.collisionTimeRes(),
                        counterColl, counterDF,
                        collIDMC,
                        posXMC, posYMC, posZMC,
                        collTimeMC,
                        isFakeCollision,
                        recoPVsPerMcColl, isPvHighestContribForMcColl);
        }

      } else {
        /// DATA
        /// Let's fill the table with all the collisions
        tableCollsBig((isRun3 ? collision.sel8() : collision.sel7()),
                      bc.runNumber(),
                      collision.posX(), collision.posY(), collision.posZ(),
                      collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(),
                      collision.numContrib(), tracksAllPerColl.size(), tracksFilteredPerColl.size(),
                      collision.chi2(),
                      bc.globalBC(),
                      bc.has_ft0() ? bc.ft0().posZ() : -999.,
                      collision.multFT0A(), collision.multFT0C(), collision.multFT0M(), collision.multFV0A(),
                      collision.collisionTime(), collision.collisionTimeRes(),
                      counterColl, counterDF,
                      collIDMC,
                      posXMC, posYMC, posZMC,
                      collTimeMC,
                      isFakeCollision,
                      recoPVsPerMcColl, isPvHighestContribForMcColl);
      }

      /// update the collision global counter
      counterColl++;
    }
  }

  /// Processing data
  void processTableDataCollsBig(CollsBigTable const& collisions,
                                soa::Filtered<TrackTableData> const& tracksFiltered,
                                TrackTableData const& tracksAll,
                                BCsWithRun3Matchings const& bcs, aod::FT0s const& ft0s)
  {

    fillCollsBigTable<false>(collisions, collisions, tracksFiltered, tracksAll, bcs, ft0s);

    /// We are at the end of the process, which is run once per DF
    /// Let's update the DF counter
    counterDF++;
  }
  PROCESS_SWITCH(QaEventTrackLiteProducer, processTableDataCollsBig, "Process data for big collision table producing", false);

  /// Processing MC
  void processTableMCCollsBig(CollsBigTableMC const& collisions,
                              McCollsWithExtra const& mcCollisions,
                              soa::Filtered<TrackTableData> const& tracksFiltered,
                              TrackTableData const& tracksAll,
                              BCsWithRun3Matchings const& bcs, aod::FT0s const& ft0s)
  {

    fillCollsBigTable<true>(collisions, mcCollisions, tracksFiltered, tracksAll, bcs, ft0s);

    /// We are at the end of the process, which is run once per DF
    /// Let's update the DF counter
    counterDF++;
  }
  PROCESS_SWITCH(QaEventTrackLiteProducer, processTableMCCollsBig, "Process MC for big collision table producing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaEventTrackLiteProducer>(cfgc)};
}
