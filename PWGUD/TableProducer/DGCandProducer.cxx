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
    registry.add("reco/TPCsignal1", "2 prong events, TPC signal of particle 1", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
    registry.add("reco/TPCsignal2", "2 prong events, TPC signal of particle 2", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
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
      LOGF(info, "DG candidate: number of PV tracks %d", collision.numContrib());
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

  // this function properly updates UDMcCollisions and UDMcParticles and returns the value
  // deltaIndex, which is needed to correct the McParticles indices
  // For a given McCollision all associated McParticles are saved
  template <typename TMcCollision, typename TMcParticles>
  void updateMcUDTables(TMcCollision const& mccol,
                        TMcParticles const& McParts,
                        int64_t& deltaIndex)
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

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMC")) {
      registry.add("mcTruth/collisions", "Number of associated collisions", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {{150, 0.0, 3.0}, {150, 0.0, 3.0}}});
    }
  }

  // process function for MC data
  // save all GRANIITTI diffractive events and the MC truth of the DG events
  void processMC(aod::McCollisions const& mccols, aod::McParticles const& mcparts,
                 UDCCs const& dgcands, UDTCs const& udtracks,
                 CCs const& collisions, BCs const& bcs, TCs const& tracks)
  {

    // loop over McCollisions and UDCCs simultaneously
    auto mccol = mccols.iteratorAt(0);
    auto dgcand = dgcands.iteratorAt(0);
    auto lastmccol = mccols.iteratorAt(mccols.size() - 1);
    auto lastdgcand = dgcands.iteratorAt(dgcands.size() - 1);

    int64_t lastSaved = -1;
    int64_t deltaIndex = 0;
    int64_t firstIndex = 0, lastIndex = 0;
    while (true) {
      // determine the next dgcand with an associated collision
      while (!dgcand.has_collision() && dgcand != lastdgcand) {
        outputMcCollsLabels(-1);
        dgcand++;
      }
      if (!dgcand.has_collision()) {
        // no dgcand left
        outputMcCollsLabels(-1);
        break;
      }

      // related mc truth
      auto dgcandCol = dgcand.collision_as<CCs>();
      if (!dgcandCol.has_mcCollision()) {
        // this collision has no MC truth
        outputMcCollsLabels(-1);
        continue;
      }
      auto mcdg = dgcandCol.mcCollision();
      auto dgmcId = mcdg.globalIndex();

      // save also all GRANIITTI diffractive events
      // keep UD Mc sorted according to AOD Mc
      auto mccolId = mccol.globalIndex();
      while (mccolId <= dgmcId) {
        auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccol.globalIndex());
        bool isGraniittiDiff = udhelpers::isGraniittiCDE(mcPartsSlice);
        bool isPythiaDiff = udhelpers::isPythiaCDE(mcPartsSlice);
        registry.get<TH1>(HIST("mcTruth/collType"))->Fill(0., 1.);
        registry.get<TH1>(HIST("mcTruth/collType"))->Fill(1., (!isPythiaDiff && !isGraniittiDiff) * 1.);
        registry.get<TH1>(HIST("mcTruth/collType"))->Fill(2., isPythiaDiff * 1.);
        registry.get<TH1>(HIST("mcTruth/collType"))->Fill(3., isGraniittiDiff * 1.);

        if (isGraniittiDiff || isPythiaDiff) {
          firstIndex = outputMcParticles.lastIndex() + 1;
          updateMcUDTables(mccol, mcPartsSlice, deltaIndex);
          lastSaved = mccolId;

          auto ivm = udhelpers::ivmGraniittiCDE(mcPartsSlice);
          registry.get<TH2>(HIST("mcTruth/IVMpt"))->Fill(ivm.M(), ivm.Perp());
        }
        if (mccol == lastmccol) {
          break;
        }
        mccol++;
        mccolId = mccol.globalIndex();
      }

      // save the MC truth of the actual dgcand
      // but check if this has not been added to the table yet
      if (lastSaved != dgmcId) {
        auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccol.globalIndex());
        firstIndex = outputMcParticles.lastIndex() + 1;
        updateMcUDTables(mccol, mcPartsSlice, deltaIndex);
        lastSaved = mccolId;

        auto ivm = udhelpers::ivmGraniittiCDE(mcPartsSlice);
        registry.get<TH2>(HIST("mcTruth/IVMpt"))->Fill(ivm.M(), ivm.Perp());
      }
      outputMcCollsLabels(outputMcCollisions.lastIndex());

      // save the mclabels of the related tracks into outputMcTrackLabels
      lastIndex = outputMcParticles.lastIndex();
      auto colTracks = udtracks.sliceByCached(aod::udtrack::udCollisionId, dgcand.globalIndex(), cache);
      for (auto colTrack : colTracks) {
        // colTrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
        auto trackId = colTrack.trackId();
        if (trackId >= 0) {
          auto track = colTrack.track_as<TCs>();
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

      // next dg candidate
      if (dgcand == lastdgcand) {
        break;
      }
      dgcand++;
    }

    // save remaining GRANIITTI diffractive events
    while (mccol != lastmccol) {
      auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccol.globalIndex());
      bool isGraniittiDiff = udhelpers::isGraniittiCDE(mcPartsSlice);
      bool isPythiaDiff = udhelpers::isPythiaCDE(mcPartsSlice);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(0., 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(1., (!isPythiaDiff && !isGraniittiDiff) * 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(2., isPythiaDiff * 1.);
      registry.get<TH1>(HIST("mcTruth/collType"))->Fill(3., isGraniittiDiff * 1.);

      if (isGraniittiDiff || isPythiaDiff) {
        updateMcUDTables(mccol, mcPartsSlice, deltaIndex);

        // update IVM versus pT
        auto ivm = udhelpers::ivmGraniittiCDE(mcPartsSlice);
        registry.get<TH2>(HIST("mcTruth/IVMpt"))->Fill(ivm.M(), ivm.Perp());
      }
      mccol++;
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
