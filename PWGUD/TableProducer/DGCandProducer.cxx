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
//
// \brief Saves relevant information of DG candidates
//
//     options:
//           DiffCuts.mNDtcoll(4)
//           DiffCuts.mMinNBCs(7)
//           DiffCuts.mMinNTracks(0)
//           DiffCuts.mMaxNTracks(10000)
//           DiffCuts.mMinNetCharge(0)
//           DiffCuts.mMaxNetCharge(0)
//           DiffCuts.mPidHypo(211)
//           DiffCuts.mMinPosz(-1000.)
//           DiffCuts.mMaxPosz(1000.)
//           DiffCuts.mMinPt(0.)
//           DiffCuts.mMaxPt(1000.)
//           DiffCuts.mMinEta(-1.)
//           DiffCuts.mMaxEta(1.)
//           DiffCuts.mMinIVM(0.)
//           DiffCuts.mMaxIVM(1000.)
//           DiffCuts.mMaxnSigmaTPC(1000.)
//           DiffCuts.mMaxnSigmaTOF(1000.)
//           DiffCutsX.mFITAmpLimits({0., 0., 0., 0., 0.})
//
//     usage: copts="--configuration json://DGCandProducerConfig.json --aod-writer-json DGCandProducerWriter.json -b"
//
//           o2-analysis-timestamp $copts |
//           o2-analysis-track-propagation $copts |
//           o2-analysis-multiplicity-table $copts |
//           o2-analysis-ft0-corrected-table $copts |
//           o2-analysis-event-selection $copts |
//           o2-analysis-trackextension $copts |
//           o2-analysis-trackselection $copts |
//           o2-analysis-pid-tpc-full $copts |
//           o2-analysis-pid-tof-base $copts |
//           o2-analysis-pid-tof-full $copts |
//           o2-analysis-ud-dgcand-producer $copts > DGCandProducer.log
//
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  20.05.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "EventFiltering/PWGUD/DGHelpers.h"
#include "PWGUD/Core/DGMCHelpers.h"
#include "PWGUD/Core/UDHelperFunctions.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandProducer {

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
  }

  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTrackCollisionIDs> outputTracksCollisionsId;

  // MC tables
  Produces<aod::UDMcCollisions> outputMcCollisions;
  Produces<aod::UDMcParticles> outputMcParticles;
  Produces<aod::UDMcTrackLabels> outputMcTrackLabels;

  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  // MC inputs
  using MCCCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MCCC = MCCCs::iterator;
  using MCTCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels, aod::TOFSignal,
                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,

                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using MCTC = MCTCs::iterator;

  // function to update UDTracks, UDTracksPID, and UDTracksExtra
  template <typename TTrack, typename TBC>
  void updateUDTrackTables(TTrack const& track, TBC const& bc)
  {
    outputTracks(track.px(), track.py(), track.pz(), track.sign(),
                 bc.globalBC(), track.trackTime(), track.trackTimeRes());
    outputTracksPID(track.tpcNSigmaEl(),
                    track.tpcNSigmaMu(),
                    track.tpcNSigmaPi(),
                    track.tpcNSigmaKa(),
                    track.tpcNSigmaPr(),
                    track.tofNSigmaEl(),
                    track.tofNSigmaMu(),
                    track.tofNSigmaPi(),
                    track.tofNSigmaKa(),
                    track.tofNSigmaPr());
    outputTracksExtra(track.itsClusterMap(),
                      track.tpcNClsFindable(),
                      track.tpcNClsFindableMinusFound(),
                      track.tpcNClsFindableMinusCrossedRows(),
                      track.tpcNClsShared(),
                      track.trdPattern(),
                      track.itsChi2NCl(),
                      track.tpcChi2NCl(),
                      track.trdChi2(),
                      track.tofChi2(),
                      track.tpcSignal(),
                      track.tofSignal(),
                      track.trdSignal(),
                      track.length(),
                      track.tofExpMom(),
                      track.detectorMap());
    outputTracksCollisionsId(outputCollisions.lastIndex());
  }

  // this function properly updates UDMcCollisions and UDMcParticles and returns the value
  // deltaIndex, which is needed to correct the McParticles indices
  // For a given McCollision all associated McParticles are saved
  template <typename TMcCollision, typename TMcParticles, typename TBC>
  void updateMcUDTables(TMcCollision const& McCol,
                        TMcParticles const& McParts,
                        TBC const& mcbc,
                        int64_t& deltaIndex)
  {
    // save McCol
    outputMcCollisions(mcbc.globalBC(),
                       McCol.generatorsID(),
                       McCol.posX(),
                       McCol.posY(),
                       McCol.posZ(),
                       McCol.t(),
                       McCol.weight(),
                       McCol.impactParameter());

    // save McParts
    // calculate conversion from old indices to new indices
    // old = mcpart.globalIndex()
    // new = old + deltaIndex
    // deltaIndex = [outputMcParticles.lastIndex() - McParts.iteratorAt(0).globalIndex() + 1]
    deltaIndex = outputMcParticles.lastIndex() - McParts.iteratorAt(0).globalIndex() + 1;
    LOGF(debug, " deltaIndex %i", deltaIndex);

    // new mother and daughter ids
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};

    // all particles of the McCollision are saved
    for (auto mcpart : McParts) {
      // correct mother and daughter IDs
      newmids.clear();
      auto oldmids = mcpart.mothersIds();
      for (uint ii = 0; ii < oldmids.size(); ii++) {
        auto newval = oldmids[ii] < 0 ? oldmids[ii] : oldmids[ii] + deltaIndex;
        LOGF(debug, " mid %i / %i", oldmids[ii], newval);
        newmids.push_back(newval);
      }
      auto olddids = mcpart.daughtersIds();
      for (uint ii = 0; ii < olddids.size(); ii++) {
        auto newval = olddids[ii] < 0 ? olddids[ii] : olddids[ii] + deltaIndex;
        LOGF(debug, " did %i / %i", olddids[ii], newval);
        newdids[ii] = newval;
      }
      LOGF(debug, " ms %i ds %i", oldmids.size(), olddids.size());

      // update UDMcParticles
      outputMcParticles(outputMcCollisions.lastIndex(),
                        mcpart.pdgCode(),
                        mcpart.statusCode(),
                        mcpart.flags(),
                        newmids,
                        newdids,
                        mcpart.weight(),
                        mcpart.px(),
                        mcpart.py(),
                        mcpart.pz(),
                        mcpart.e());
    }
  }

  // process function for real data
  void processData(CC const& collision,
                   BCs const& bcs,
                   TCs& tracks,
                   FWs& fwdtracks,
                   aod::Zdcs& zdcs,
                   aod::FT0s& ft0s,
                   aod::FV0As& fv0as,
                   aod::FDDs& fdds)
  {
    // nominal BC
    auto bc = collision.bc_as<BCs>();

    // obtain slice of compatible BCs
    auto bcRange = compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, tracks, fwdtracks);

    // save DG candidates
    if (isDGEvent == 0) {
      LOGF(info, "  Data: good collision!");

      // update DG candidates tables
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       collision.numContrib(), netCharge(tracks),
                       rPVtrwTOF(tracks, collision.numContrib()),
                       0., 0., 0., 0., 0);

      // update DGTracks tables
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          updateUDTrackTables(track, bc);
        }
      }
    }
  }

  PROCESS_SWITCH(DGCandProducer, processData, "Process real data", false);

  Preslice<MCTCs> tracksPerCollision = aod::track::collisionId;
  Preslice<FWs> fwdTracksPerCollision = aod::fwdtrack::collisionId;

  // process function for MC data
  void processMc(aod::McCollision const& McCol,
                 aod::McParticles const& McParts,
                 MCCCs const& collisions,
                 BCs const& bcs,
                 MCTCs const& tracks,
                 FWs const& fwdtracks,
                 aod::Zdcs const& zdcs,
                 aod::FT0s const& ft0s,
                 aod::FV0As const& fv0as,
                 aod::FDDs const& fdds) //)
  {
    for (auto McPart : McParts) {
      LOGF(debug, "McCol %i McPart %i", McCol.globalIndex(), McPart.globalIndex());
    }

    // is this a central diffractive event?
    // by default it is assumed to be a MB event
    bool isPythiaDiff = isPythiaCDE(McParts);
    bool isGraniittiDiff = isGraniittiCDE(McParts);
    LOGF(debug, "mcCol %i type %i / %i / %i", (int)McCol.globalIndex(), !isPythiaDiff && !isGraniittiDiff, isPythiaDiff, isGraniittiDiff);

    // MC BC
    auto mcbc = McCol.bc_as<BCs>();

    // save MCTruth of all diffractive events
    bool mcColIsSaved = false;
    int64_t deltaIndex = 0;
    auto nMcParts0 = outputMcParticles.lastIndex();
    if (isPythiaDiff || isGraniittiDiff) {
      // update tables UDMcCollisions and UDMcParticles
      updateMcUDTables(McCol, McParts, mcbc, deltaIndex);
      mcColIsSaved = true;
    }

    // loop over all reconstructed collisions associated with McCol
    for (auto collision : collisions) {
      // get the tracks belonging to collision
      auto collisionTracks = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      auto collisionFwdTracks = fwdtracks.sliceBy(fwdTracksPerCollision, collision.globalIndex());
      LOGF(debug, "  tracks %i / %i", (int)collisionTracks.size(), collisionFwdTracks.size());

      auto bc = collision.bc_as<BCs>();
      LOGF(debug, "  BC mc %i reco %i", mcbc.globalBC(), bc.globalBC());

      // is this a collision to be saved?
      // obtain slice of compatible BCs
      auto bcRange = MCcompatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());

      // apply DG selection
      auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, collisionTracks, collisionFwdTracks);
      LOGF(debug, "  isDG %i", (int)isDGEvent);

      // save information of DG events
      if (isDGEvent == 0) {
        LOGF(info, "  MC: good collision!");

        // update UDMcCollisions and UDMcParticles if not already done
        if (!mcColIsSaved) {
          updateMcUDTables(McCol, McParts, mcbc, deltaIndex);
          mcColIsSaved = true;
        }

        // UDCollisions
        outputCollisions(bc.globalBC(), bc.runNumber(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         collision.numContrib(), netCharge(tracks),
                         rPVtrwTOF(collisionTracks, collision.numContrib()),
                         0., 0., 0., 0., 0);

        // UDTracks, UDTrackCollisionID, UDTracksExtras, UDMcTrackLabels
        for (auto& track : collisionTracks) {
          // but save only the Primary Vertex tracks
          if (track.isPVContributor()) {
            updateUDTrackTables(track, bc);

            // properly correct the index into the UDMcParticles tables with deltaIndex
            auto newval = track.mcParticleId() < 0 ? track.mcParticleId() : track.mcParticleId() + deltaIndex;
            // only associations with McParticles belonging to the actual McCollision are supported
            if ((newval < nMcParts0) || (newval > outputMcParticles.lastIndex())) {
              LOGF(info, "<ATTENTION> UDMcParticles index out of range %i (%i - %i)", newval, nMcParts0 + 1, outputMcParticles.lastIndex());
              newval = -1;
            }
            outputMcTrackLabels(newval,
                                track.mcMask());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DGCandProducer, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandProducer>(cfgc, TaskName{"dgcandproducer"}),
  };
}
