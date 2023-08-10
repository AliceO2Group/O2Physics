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
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  20.05.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/DGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandProducer {
  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDFwdTracks> outputFwdTracks;
  Produces<aod::UDFwdTracksExtra> outputFwdTracksExtra;

  // MC tables
  Produces<aod::UDMcCollisions> outputMcCollisions;
  Produces<aod::UDMcParticles> outputMcParticles;
  Produces<aod::UDMcTrackLabels> outputMcTrackLabels;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFbeta,
                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  // MC inputs
  using MCCCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MCCC = MCCCs::iterator;
  using MCTCs = soa::Join<aod::Tracks, aod::TracksExtra, /*aod::TracksCov,*/ aod::TracksDCA, aod::TrackSelection,
                          aod::McTrackLabels,
                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                          aod::TOFSignal, aod::pidTOFbeta,
                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using MCTC = MCTCs::iterator;

  // function to update UDFwdTracks, UDFwdTracksExtra
  template <typename TFwdTrack>
  void updateUDFwdTrackTables(TFwdTrack const& fwdtrack, uint64_t const& bcnum)
  {
    outputFwdTracks(outputCollisions.lastIndex(),
                    fwdtrack.px(), fwdtrack.py(), fwdtrack.pz(), fwdtrack.sign(),
                    bcnum, fwdtrack.trackTime(), fwdtrack.trackTimeRes());
    outputFwdTracksExtra(fwdtrack.nClusters(),
                         fwdtrack.pDca(),
                         fwdtrack.rAtAbsorberEnd(),
                         fwdtrack.chi2(),
                         fwdtrack.chi2MatchMCHMID(),
                         fwdtrack.mchBitMap(),
                         fwdtrack.midBitMap(),
                         fwdtrack.midBoards());
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
    LOGF(debug, "  deltaIndex (%d) = lastIndex (%d) - McPartsfirst (%d) + 1", deltaIndex, outputMcParticles.lastIndex(), McParts.iteratorAt(0).globalIndex());

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

  SliceCache cache;
  PresliceUnsorted<aod::McParticles> mcPartsPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<MCCCs> collisionsPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  Preslice<MCTCs> tracksPerCollision = aod::track::collisionId;
  Preslice<FWs> fwdTracksPerCollision = aod::fwdtrack::collisionId;

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack>
  void updateUDTrackTables(int64_t lastIndex, TTrack const& track, uint64_t const& bcnum)
  {
    outputTracks(lastIndex,
                 track.px(), track.py(), track.pz(), track.sign(),
                 bcnum, track.trackTime(), track.trackTimeRes());

    // float sigmaY = track.sigmaY();
    // float sigmaZ = track.sigmaZ();
    float sigmaY = -1.;
    float sigmaZ = -1.;
    outputTracksCov(track.x(), track.y(), track.z(), sigmaY, sigmaZ);

    outputTracksDCA(track.dcaZ(), track.dcaXY());
    outputTracksPID(track.tpcNSigmaEl(),
                    track.tpcNSigmaMu(),
                    track.tpcNSigmaPi(),
                    track.tpcNSigmaKa(),
                    track.tpcNSigmaPr(),
                    track.beta(),
                    track.betaerror(),
                    track.tofNSigmaEl(),
                    track.tofNSigmaMu(),
                    track.tofNSigmaPi(),
                    track.tofNSigmaKa(),
                    track.tofNSigmaPr());
    outputTracksExtra(track.tpcInnerParam(),
                      track.itsClusterMap(),
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
    outputTracksFlag(track.has_collision(),
                     track.isPVContributor());
  }

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processData")) {
      registry.add("data/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}});
      registry.add("data/pt1Vspt2", "2 prong events, p_{T} versus p_{T}", {HistType::kTH2F, {{100, -3., 3.}, {100, -3., 3.0}}});
      registry.add("data/TPCsignal1", "2 prong events, TPC signal of particle 1", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
      registry.add("data/TPCsignal2", "2 prong events, TPC signal of particle 2", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
      registry.add("data/sig1VsSig2TPC", "2 prong events, TPC signal versus TPC signal", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}});
    }

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("mcTruth/collisions", "Number of associated collisions", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {{150, 0.0, 3.0}, {150, 0.0, 3.0}}});
    }
  }

  // process function for real data
  void processData(CC const& collision, BCs const& bcs, TCs& tracks, FWs& fwdtracks,
                   aod::Zdcs& zdcs, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    LOGF(debug, "<DGCandProducer>  collision %d", collision.globalIndex());
    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    auto bc = collision.foundBC_as<BCs>();
    LOGF(debug, "<DGCandProducer>  BC id %d", bc.globalBC());

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
    LOGF(debug, "<DGCandProducer>  Size of bcRange %d", bcRange.size());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, tracks, fwdtracks);

    // save DG candidates
    registry.get<TH1>(HIST("data/Stat"))->Fill(0., 1.);
    registry.get<TH1>(HIST("data/Stat"))->Fill(isDGEvent + 1, 1.);
    if (isDGEvent == 0) {
      LOGF(debug, "<DGCandProducer>  Data: good collision!");

      // fill FITInfo
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bc.globalBC(), bcs, ft0s, fv0as, fdds);

      // update DG candidates tables
      auto rtrwTOF = udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib());
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       collision.numContrib(), udhelpers::netCharge<true>(tracks),
                       rtrwTOF);
      outputCollisionsSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C,
                           fitInfo.triggerMaskFT0,
                           fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC,
                           fitInfo.triggerMaskFDD,
                           fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                           fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                           fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                           fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);

      // update DGTracks tables
      for (auto& track : tracks) {
        updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
      }

      // update DGFwdTracks tables
      for (auto& fwdtrack : fwdtracks) {
        updateUDFwdTrackTables(fwdtrack, bc.globalBC());
      }

      // fill UDZdcs
      if (bc.has_zdc()) {
        auto zdc = bc.zdc();
        auto enes = std::vector(zdc.energy().begin(), zdc.energy().end());
        auto chEs = std::vector(zdc.channelE().begin(), zdc.channelE().end());
        auto amps = std::vector(zdc.amplitude().begin(), zdc.amplitude().end());
        auto times = std::vector(zdc.time().begin(), zdc.time().end());
        auto chTs = std::vector(zdc.channelT().begin(), zdc.channelT().end());
        outputZdcs(outputCollisions.lastIndex(), enes, chEs, amps, times, chTs);
      }

      // produce TPC signal histograms for 2-track events
      if (collision.numContrib() == 2) {
        auto cnt = 0;
        float pt1 = 0., pt2 = 0.;
        float signalTPC1 = 0., signalTPC2 = 0.;
        for (auto tr : tracks) {
          if (tr.isPVContributor()) {
            cnt++;
            switch (cnt) {
              case 1:
                pt1 = tr.pt() * tr.sign();
                signalTPC1 = tr.tpcSignal();
                break;
              case 2:
                pt2 = tr.pt() * tr.sign();
                signalTPC2 = tr.tpcSignal();
            }
            LOGF(debug, "<DGCandProducer>    track[%d] %d pT %f ITS %d TPC %d TRD %d TOF %d",
                 cnt, tr.isGlobalTrack(), tr.pt(), tr.itsNCls(), tr.tpcNClsCrossedRows(), tr.hasTRD(), tr.hasTOF());
          }
        }
        registry.get<TH2>(HIST("data/pt1Vspt2"))->Fill(pt1, pt2);
        registry.get<TH2>(HIST("data/TPCsignal1"))->Fill(pt1, signalTPC1);
        registry.get<TH2>(HIST("data/TPCsignal2"))->Fill(pt2, signalTPC2);
        registry.get<TH2>(HIST("data/sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);
      }
    }
  }
  PROCESS_SWITCH(DGCandProducer, processData, "Process real data", false);

  // process function for MC data
  void processMc(aod::McCollisions const& McCols,
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
    LOGF(info, "Number of McCollisions %d", McCols.size());
    for (auto McCol : McCols) {
      // MC BC and Particles
      auto mcbc = McCol.bc_as<BCs>();
      auto mcPartsSlice = McParts.sliceBy(mcPartsPerMcCollision, McCol.globalIndex());

      // is this a central diffractive event?
      // by default it is assumed to be a MB event
      bool isPythiaDiff = udhelpers::isPythiaCDE(mcPartsSlice);
      bool isGraniittiDiff = udhelpers::isGraniittiCDE(mcPartsSlice);

      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(0., 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(1., (!isPythiaDiff && !isGraniittiDiff) * 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(2., isPythiaDiff * 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(3., isGraniittiDiff * 1.);

      // save MCTruth of all diffractive events
      bool mcColIsSaved = false;
      int64_t deltaIndex = 0;
      auto nMcParts0 = outputMcParticles.lastIndex();
      if (isPythiaDiff || isGraniittiDiff) {
        // update tables UDMcCollisions and UDMcParticles
        updateMcUDTables(McCol, mcPartsSlice, mcbc, deltaIndex);
        mcColIsSaved = true;

        // update IVM versus pT
        auto ivm = udhelpers::ivmGraniittiCDE(mcPartsSlice);
        registry.get<TH2>(HIST("mcTruth/IVMpt"))->Fill(ivm.M(), ivm.Perp());
      }

      // loop over all reconstructed collisions associated with McCol
      auto collisionsSlice = collisions.sliceBy(collisionsPerMcCollision, McCol.globalIndex());
      registry.get<TH1>(HIST("mcTruth/collisions"))->Fill(collisionsSlice.size(), 1.);

      for (auto collision : collisionsSlice) {
        // get the tracks belonging to collision
        auto collisionTracks = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        auto collisionFwdTracks = fwdtracks.sliceByCached(aod::fwdtrack::collisionId, collision.globalIndex(), cache);

        auto bc = collision.bc_as<BCs>();

        // is this a collision to be saved?
        // obtain slice of compatible BCs
        auto bcRange = udhelpers::MCcompatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());

        // apply DG selection
        auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, collisionTracks, collisionFwdTracks);

        // save information of DG events
        if (isDGEvent == 0) {

          // update UDMcCollisions and UDMcParticles if not already done
          if (!mcColIsSaved) {
            updateMcUDTables(McCol, mcPartsSlice, mcbc, deltaIndex);
            mcColIsSaved = true;
          }

          // UDCollisions
          auto rtrwTOF = udhelpers::rPVtrwTOF<true>(collisionTracks, collision.numContrib());
          outputCollisions(bc.globalBC(), bc.runNumber(),
                           collision.posX(), collision.posY(), collision.posZ(),
                           collision.numContrib(), udhelpers::netCharge<true>(tracks),
                           rtrwTOF);

          // UDTracks, UDTrackCollisionID, UDTracksExtras, UDMcTrackLabels
          for (auto& track : collisionTracks) {
            updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());

            // properly correct the index into the UDMcParticles tables with deltaIndex
            auto newval = track.mcParticleId() < 0 ? track.mcParticleId() : track.mcParticleId() + deltaIndex;
            // only associations with McParticles belonging to the actual McCollision are supported
            if ((newval < nMcParts0) || (newval > outputMcParticles.lastIndex())) {
              LOGF(debug, "  <ATTENTION> UDMcParticles index out of range %i (%i - %i)", newval, nMcParts0 + 1, outputMcParticles.lastIndex());
              newval = -1;
            }
            outputMcTrackLabels(newval, track.mcMask());
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
