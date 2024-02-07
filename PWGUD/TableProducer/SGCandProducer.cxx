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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/SGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SGCandProducer {
  // get an SGCutparHolder
  SGCutParHolder sameCuts = SGCutParHolder(); // SGCutparHolder
  Configurable<SGCutParHolder> SGCuts{"SGCuts", {}, "SG event cuts"};

  // SG selector
  SGSelector sgSelector;

  // data tables
  Produces<aod::SGCollisions> outputSGCollisions;
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
                      track.itsClusterSizes(),
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
    sameCuts = (SGCutParHolder)SGCuts;

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
    LOGF(debug, "<SGCandProducer>  collision %d", collision.globalIndex());
    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    auto bc = collision.foundBC_as<BCs>();
    LOGF(debug, "<SGCandProducer>  BC id %d", bc.globalBC());

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, sameCuts.NDtcoll(), bcs, sameCuts.minNBCs());
    LOGF(debug, "<SGCandProducer>  Size of bcRange %d", bcRange.size());

    // apply SG selection
    //    auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, tracks, fwdtracks);
    auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, tracks, fwdtracks);

    // Check if it's a SingleGap event for sideC
    registry.get<TH1>(HIST("reco/Stat"))->Fill(0., 1.);
    registry.get<TH1>(HIST("reco/Stat"))->Fill(isSGEvent + 1, 1.);
    if (isSGEvent <= 2) {
      LOGF(debug, "<SGCandProducer>  Data: good collision!");

      // fill FITInfo
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bc.globalBC(), bcs, ft0s, fv0as, fdds);

      // update SG candidates tables
      auto rtrwTOF = udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib());
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       collision.numContrib(), udhelpers::netCharge<true>(tracks),
                       rtrwTOF);
      outputSGCollisions(isSGEvent);
      outputCollisionsSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C,
                           fitInfo.triggerMaskFT0,
                           fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC,
                           fitInfo.triggerMaskFDD,
                           fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                           fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                           fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                           fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      outputCollsLabels(collision.globalIndex());

      // update SGTracks tables
      for (auto& track : tracks) {
        updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
      }

      // update SGFwdTracks tables
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
      LOGF(debug, "SG candidate: number of PV tracks %d", collision.numContrib());
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
            LOGF(debug, "<SGCandProducer>    track[%d] %d pT %f ITS %d TPC %d TRD %d TOF %d",
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

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<SGCandProducer>(cfgc, TaskName{"sgcandproducer"})};

  return workflow;
}
