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

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
  }

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

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"pt1Vspt2", "#pt1Vspt2", {HistType::kTH2F, {{100, -3., 3.}, {100, -3., 3.0}}}},
      {"TPCsignal1", "#TPCsignal1", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}}},
      {"TPCsignal2", "#TPCsignal2", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}}},
      {"sig1VsSig2TPC", "#sig1VsSig2TPC", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}}},
    }};

  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  // MC inputs
  using MCCCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MCCC = MCCCs::iterator;
  using MCTCs = soa::Join<aod::Tracks, aod::TracksExtra, /*aod::TracksCov,*/ aod::TracksDCA, aod::TrackSelection,
                          aod::McTrackLabels,
                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                          aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using MCTC = MCTCs::iterator;

  // extract FIT information
  void getFITinfo(upchelpers::FITInfo& info, uint64_t const& bcnum, BCs const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // find bc with globalBC = bcnum
    Partition<BCs> selbc = aod::bc::globalBC == bcnum;
    selbc.bindTable(bcs);

    // if BC exists then update FIT information for this BC
    if (selbc.size() > 0) {
      auto bc = selbc.begin();

      // FT0
      if (bc.has_foundFT0()) {
        auto ft0 = ft0s.iteratorAt(bc.foundFT0Id());
        info.timeFT0A = ft0.timeA();
        info.timeFT0C = ft0.timeC();
        const auto& ampsA = ft0.amplitudeA();
        const auto& ampsC = ft0.amplitudeC();
        info.ampFT0A = 0.;
        for (auto amp : ampsA) {
          info.ampFT0A += amp;
        }
        info.ampFT0C = 0.;
        for (auto amp : ampsC) {
          info.ampFT0C += amp;
        }
        info.triggerMaskFT0 = ft0.triggerMask();
      }

      // FV0A
      if (bc.has_foundFV0()) {
        auto fv0a = fv0as.iteratorAt(bc.foundFV0Id());
        info.timeFV0A = fv0a.time();
        const auto& amps = fv0a.amplitude();
        info.ampFV0A = 0.;
        for (auto amp : amps) {
          info.ampFV0A += amp;
        }
        info.triggerMaskFV0A = fv0a.triggerMask();
      }

      // FDD
      if (bc.has_foundFDD()) {
        auto fdd = fdds.iteratorAt(bc.foundFDDId());
        info.timeFDDA = fdd.timeA();
        info.timeFDDC = fdd.timeC();
        const auto& ampsA = fdd.chargeA();
        const auto& ampsC = fdd.chargeC();
        info.ampFDDA = 0.;
        for (auto amp : ampsA) {
          info.ampFDDA += amp;
        }
        info.ampFDDC = 0.;
        for (auto amp : ampsC) {
          info.ampFDDC += amp;
        }
        info.triggerMaskFDD = fdd.triggerMask();
      }
    }

    // fill BG and BB flags in adjacent BCs [-16, 15]
    // compute range to check
    auto minbc = bcnum - 16;
    auto maxbc = bcnum + 15;
    Partition<BCs> bcrange = aod::bc::globalBC >= minbc && aod::bc::globalBC <= maxbc;
    bcrange.bindTable(bcs);

    // loop over bcrange and check
    for (auto const& bc2u : bcrange) {

      // 0 <= bit <= 31
      auto bit = bc2u.globalBC() - minbc;
      if (!bc2u.selection_bit(evsel::kNoBGT0A))
        SETBIT(info.BGFT0Apf, bit);
      if (!bc2u.selection_bit(evsel::kNoBGT0C))
        SETBIT(info.BGFT0Cpf, bit);
      if (bc2u.selection_bit(evsel::kIsBBT0A))
        SETBIT(info.BBFT0Apf, bit);
      if (bc2u.selection_bit(evsel::kIsBBT0C))
        SETBIT(info.BBFT0Cpf, bit);
      if (!bc2u.selection_bit(evsel::kNoBGV0A))
        SETBIT(info.BGFV0Apf, bit);
      if (bc2u.selection_bit(evsel::kIsBBV0A))
        SETBIT(info.BBFV0Apf, bit);
      if (!bc2u.selection_bit(evsel::kNoBGFDA))
        SETBIT(info.BGFDDApf, bit);
      if (!bc2u.selection_bit(evsel::kNoBGFDC))
        SETBIT(info.BGFDDCpf, bit);
      if (bc2u.selection_bit(evsel::kIsBBFDA))
        SETBIT(info.BBFDDApf, bit);
      if (bc2u.selection_bit(evsel::kIsBBFDC))
        SETBIT(info.BBFDDCpf, bit);
    }
  }

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack>
  void updateUDTrackTables(TTrack const& track, uint64_t const& bcnum)
  {
    outputTracks(outputCollisions.lastIndex(),
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
    if (isDGEvent == 0) {
      LOGF(debug, "<DGCandProducer>  Data: good collision!");

      // fill FITInfo
      upchelpers::FITInfo fitInfo{};
      getFITinfo(fitInfo, bc.globalBC(), bcs, ft0s, fv0as, fdds);

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

      // update DGTracks tables
      for (auto& track : tracks) {
        updateUDTrackTables(track, bc.globalBC());
      }

      // update DGFwdTracks tables
      for (auto& fwdtrack : fwdtracks) {
        updateUDFwdTrackTables(fwdtrack, bc.globalBC());
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
        registry.get<TH2>(HIST("pt1Vspt2"))->Fill(pt1, pt2);
        registry.get<TH2>(HIST("TPCsignal1"))->Fill(pt1, signalTPC1);
        registry.get<TH2>(HIST("TPCsignal2"))->Fill(pt2, signalTPC2);
        registry.get<TH2>(HIST("sig1VsSig2TPC"))->Fill(signalTPC1, signalTPC2);
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
    bool isPythiaDiff = udhelpers::isPythiaCDE(McParts);
    bool isGraniittiDiff = udhelpers::isGraniittiCDE(McParts);
    LOGF(debug, "mcCol %i type %i / %i / %i", (int)McCol.globalIndex(), !isPythiaDiff && !isGraniittiDiff, isPythiaDiff, isGraniittiDiff);
    /*
    // mctruth
    int mctruth = -1;
    if (isPythiaDiff) {
      mctruth = 1;
    } else if (isGraniittiDiff) {
      mctruth = 2;
    }
    */

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
      auto bcRange = udhelpers::MCcompatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());

      // apply DG selection
      auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, collisionTracks, collisionFwdTracks);
      LOGF(debug, "  isDG %i", (int)isDGEvent);

      // save information of DG events
      if (isDGEvent == 0) {
        LOGF(debug, "  MC: good collision!");

        // update UDMcCollisions and UDMcParticles if not already done
        if (!mcColIsSaved) {
          updateMcUDTables(McCol, McParts, mcbc, deltaIndex);
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
          // but save only the Primary Vertex tracks
          updateUDTrackTables(track, bc.globalBC());

          // properly correct the index into the UDMcParticles tables with deltaIndex
          auto newval = track.mcParticleId() < 0 ? track.mcParticleId() : track.mcParticleId() + deltaIndex;
          // only associations with McParticles belonging to the actual McCollision are supported
          if ((newval < nMcParts0) || (newval > outputMcParticles.lastIndex())) {
            LOGF(info, "<ATTENTION> UDMcParticles index out of range %i (%i - %i)", newval, nMcParts0 + 1, outputMcParticles.lastIndex());
            newval = -1;
          }
          outputMcTrackLabels(newval, track.mcMask());
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
