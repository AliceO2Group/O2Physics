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
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDFwdTracks> outputFwdTracks;
  Produces<aod::UDFwdTracksExtra> outputFwdTracksExtra;
  Produces<aod::UDTracksLabels> outputTracksLabel;

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
    outputTracksLabel(track.globalIndex());
  }

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    // add histograms for the different process functions
    registry.add("reco/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}});
    registry.add("reco/pt1Vspt2", "2 prong events, p_{T} versus p_{T}", {HistType::kTH2F, {{100, -3., 3.}, {100, -3., 3.0}}});
    registry.add("reco/TPCsignal1", "2 prong events, TPC signal versus p_{T} of particle 1", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
    registry.add("reco/TPCsignal2", "2 prong events, TPC signal versus p_{T} of particle 2", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
    registry.add("reco/sig1VsSig2TPC", "2 prong events, TPC signal versus TPC signal", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}});
  }

  // process function for real data
  void process(CC const& collision, BCs const& bcs, TCs& tracks, FWs& fwdtracks,
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
    registry.get<TH1>(HIST("reco/Stat"))->Fill(0., 1.);
    registry.get<TH1>(HIST("reco/Stat"))->Fill(isDGEvent + 1, 1.);
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
      outputCollsLabels(collision.globalIndex());

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
      LOGF(debug, "DG candidate: number of PV tracks %d", collision.numContrib());
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
        registry.get<TH2>(HIST("reco/pt1Vspt2"))->Fill(pt1, pt2);
        registry.get<TH2>(HIST("reco/TPCsignal1"))->Fill(pt1, signalTPC1);
        registry.get<TH2>(HIST("reco/TPCsignal2"))->Fill(pt2, signalTPC2);
        registry.get<TH2>(HIST("reco/sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);
      }
    }
  }
};

struct McDGCandProducer {
  // MC tables
  Produces<aod::UDMcCollisions> outputMcCollisions;
  Produces<aod::UDMcParticles> outputMcParticles;
  Produces<aod::UDMcCollsLabels> outputMcCollsLabels;
  Produces<aod::UDMcTrackLabels> outputMcTrackLabels;

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using UDCCs = soa::Join<aod::UDCollisions, aod::UDCollsLabels>;
  using UDTCs = soa::Join<aod::UDTracks, aod::UDTracksLabels>;

  // prepare slices
  SliceCache cache;
  PresliceUnsorted<aod::McParticles> mcPartsPerMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<UDTCs> udtracksPerUDCollision = aod::udtrack::udCollisionId;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  template <typename TMcCollision>
  void updateUDMcCollisions(TMcCollision const& mccol)
  {
    // save mccol
    outputMcCollisions(mccol.bcId(),
                       mccol.generatorsID(),
                       mccol.posX(),
                       mccol.posY(),
                       mccol.posZ(),
                       mccol.t(),
                       mccol.weight(),
                       mccol.impactParameter());
  }

  template <typename TMcParticle>
  void updateUDMcParticle(TMcParticle const& McPart, int64_t& deltaIndex, int64_t McCollisionId)
  {
    // save McPart
    // mother and daughter indices are set to -1
    // ATTENTION: this can be improved to also include mother and daughter indices
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};

    // update UDMcParticles
    outputMcParticles(McCollisionId,
                      McPart.pdgCode(),
                      McPart.statusCode(),
                      McPart.flags(),
                      newmids,
                      newdids,
                      McPart.weight(),
                      McPart.px(),
                      McPart.py(),
                      McPart.pz(),
                      McPart.e());
  }

  template <typename TMcParticles>
  void updateUDMcParticles(TMcParticles const& McParts, int64_t& deltaIndex, int64_t McCollisionId)
  {
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
      outputMcParticles(McCollisionId,
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

  template <typename TTrack>
  void updateUDMcTrackLabel(TTrack const& udtrack, int64_t const& deltaIndex)
  {
    // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
    auto trackId = udtrack.trackId();
    if (trackId >= 0) {
      auto track = udtrack.template track_as<TCs>();
      auto mcTrackId = track.mcParticleId();
      if (mcTrackId >= 0) {
        auto udMcTrackId = mcTrackId + deltaIndex;
        outputMcTrackLabels(udMcTrackId, track.mcMask());
      } else {
        outputMcTrackLabels(-1, track.mcMask());
      }
    } else {
      outputMcTrackLabels(-1, -1);
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabels(TTrack const& udtracks, int64_t const& deltaIndex)
  {
    // loop over all tracks
    for (auto udtrack : udtracks) {
      // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
      auto trackId = udtrack.trackId();
      if (trackId >= 0) {
        auto track = udtrack.template track_as<TCs>();
        auto mcTrackId = track.mcParticleId();
        if (mcTrackId >= 0) {
          auto udMcTrackId = mcTrackId + deltaIndex;
          outputMcTrackLabels(udMcTrackId, track.mcMask());
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, -1);
      }
    }
  }

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMC")) {
      registry.add("mcTruth/collisions", "Number of associated collisions", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {{150, 0.0, 3.0}, {150, 0.0, 3.0}}});
    }
  }

  // process function for MC data
  // save the MC truth of all events of interest and of the DG events
  void processMC(aod::McCollisions const& mccols, aod::McParticles const& mcparts,
                 UDCCs const& dgcands, UDTCs const& udtracks,
                 CCs const& collisions, BCs const& bcs, TCs const& tracks)
  {
    LOGF(info, "Number of McCollisions %d", mccols.size());
    LOGF(info, "Number of DG candidates %d", dgcands.size());
    LOGF(info, "Number of UD tracks %d", udtracks.size());

    // use a hash table to keep track of the McCollisions which have been added to the UDMcCollision table
    // {McCollisionId : udMcCollisionId}
    std::map<int64_t, int64_t> mcColIsSaved;
    /*
      example code for a list with dgcnad and mccol

        std::vector<std::vector<int>> qq;
        qq.push_back(std::vector<int>{1,1});
        qq.push_back(std::vector<int>{2,2});
        qq.push_back(std::vector<int>{3,-1});
        qq.push_back(std::vector<int>{5,-1});
        qq.push_back(std::vector<int>{7,9});
        qq.push_back(std::vector<int>{-1,3});
        qq.push_back(std::vector<int>{-1,4});
        qq.push_back(std::vector<int>{-1,5});

        std::sort(qq.begin(), qq.end(),
          [](const std::vector<int>& a, const std::vector<int>& b) {
            if (a[1] != -1 && b[1] != -1) {
              return a[1] < b[1];
            } else {
              return a[0] < b[0];
            }
        });

        for (auto q : qq) {
          LOGF(info, " q[0] %d q [1] %d", q[0], q[1]);
        }
    */

    // loop over McCollisions and UDCCs simultaneously
    auto mccol = mccols.iteratorAt(0);
    auto dgcand = dgcands.iteratorAt(0);
    auto lastmccol = mccols.iteratorAt(mccols.size() - 1);
    auto lastdgcand = dgcands.iteratorAt(dgcands.size() - 1);

    int64_t deltaIndex = 0;

    // advance dgcand and mccol until both are AtEnd
    int64_t mccolId = mccol.globalIndex();
    int64_t mcdgId = -1;
    auto dgcandAtEnd = dgcand == lastdgcand;
    auto mccolAtEnd = mccol == lastmccol;
    bool goon = true;
    while (goon) {
      // check if dgcand has an associated McCollision
      if (!dgcand.has_collision()) {
        mcdgId = -1;
      } else {
        auto dgcandCol = dgcand.collision_as<CCs>();
        if (!dgcandCol.has_mcCollision()) {
          mcdgId = -1;
        } else {
          mcdgId = dgcandCol.mcCollision().globalIndex();
        }
      }
      LOGF(info, "\nstart of loop mcdgId %d mccolId %d", mcdgId, mccolId);

      // two cases to consider
      // 1. the event to process is a dgcand. In this case the Mc tables as well as the McLabel tables are updated
      // 2. the event to process is an event of interest. In this case only the Mc tables are updated
      if ((!dgcandAtEnd && !mccolAtEnd && (mcdgId <= mccolId)) || mccolAtEnd) {
        // this is case 1.
        LOGF(info, "doing case 1 with mcdgId %d", mcdgId);

        // update UDMcCollisions and UDMcColsLabels (for each UDCollision -> UDMcCollisions)
        // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
        // get dgcand tracks
        auto dgTracks = udtracks.sliceByCached(aod::udtrack::udCollisionId, dgcand.globalIndex(), cache);

        // If the dgcand has an associated McCollision then the McCollision and all associated
        // McParticles are saved
        if (mcdgId >= 0) {
          LOGF(info, "  saving McCollision %d", mcdgId);
          if (mcColIsSaved.find(mcdgId) == mcColIsSaved.end()) {
            // update UDMcCollisions
            auto dgcandMcCol = dgcand.collision_as<CCs>().mcCollision();
            updateUDMcCollisions(dgcandMcCol);
            mcColIsSaved[mcdgId] = outputMcCollisions.lastIndex();
          }

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(mcColIsSaved[mcdgId]);

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mcdgId);
          updateUDMcParticles(mcPartsSlice, deltaIndex, mcColIsSaved[mcdgId]);

          // update UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          updateUDMcTrackLabels(dgTracks, deltaIndex);

        } else {
          // If the dgcand has no associated McCollision then only the McParticles which are associated
          // with the tracks of the dgcand are saved
          LOGF(info, "  saving McCollision %d", -1);

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(-1);

          // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          // loop over tracks of dgcand
          for (auto dgtrack : dgTracks) {
            if (dgtrack.has_track()) {
              auto track = dgtrack.track_as<TCs>();
              if (track.has_mcParticle()) {
                auto mcPart = track.mcParticle();
                updateUDMcParticle(mcPart, deltaIndex, -1);
                updateUDMcTrackLabel(dgtrack, deltaIndex);
              } else {
                outputMcTrackLabels(-1, track.mcMask());
              }
            } else {
              outputMcTrackLabels(-1, -1);
            }
          }
        }
        // advance dgcand
        if (dgcand != lastdgcand) {
          dgcand++;
        } else {
          dgcandAtEnd = true;
        }
      } else {
        // this is case 2.
        LOGF(info, "doing case 2");

        // update UDMcCollisions and UDMcParticles
        if (mcColIsSaved.find(mccolId) == mcColIsSaved.end()) {
          LOGF(info, "  saving McCollision %d", mccolId);
          // update UDMcCollisions
          updateUDMcCollisions(mccol);
          mcColIsSaved[mccolId] = outputMcCollisions.lastIndex();

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccolId);
          updateUDMcParticles(mcPartsSlice, deltaIndex, mcColIsSaved[mccolId]);
        }

        // advance mccol
        if (mccol != lastmccol) {
          mccol++;
          mccolId = mccol.globalIndex();
        } else {
          mccolAtEnd = true;
        }
      }

      goon = !dgcandAtEnd || !mccolAtEnd;
      LOGF(info, "end of loop mcdgId %d mccolId %d", mcdgId, mccolId);
    }
  }
  PROCESS_SWITCH(McDGCandProducer, processMC, "Produce MC tables", false);

  void processDummy(aod::Collisions const& collisions)
  {
    // do nothing
    LOGF(info, "Running dummy process function!");
  }
  PROCESS_SWITCH(McDGCandProducer, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DGCandProducer>(cfgc, TaskName{"dgcandproducer"}),
    adaptAnalysisTask<McDGCandProducer>(cfgc, TaskName{"mcdgcandproducer"})};

  return workflow;
}
