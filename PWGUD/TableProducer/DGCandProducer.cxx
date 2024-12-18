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

#include <vector>
#include <string>
#include <map>
#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/DGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandProducer {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher mRateFetcher;
  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  Configurable<bool> saveAllTracks{"saveAllTracks", true, "save only PV contributors or all tracks associated to a collision"};
  Configurable<bool> fillFIThistos{"fillFIThistos", false, "fill the histograms with the FIT amplitudes"};

  // DG selector
  DGSelector dgSelector;

  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDCollisionSelExtras> outputCollisionSelExtras;
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDZdcsReduced> outputZdcsReduced;
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
    outputFwdTracksExtra(fwdtrack.trackType(),
                         fwdtrack.nClusters(),
                         fwdtrack.pDca(),
                         fwdtrack.rAtAbsorberEnd(),
                         fwdtrack.chi2(),
                         fwdtrack.chi2MatchMCHMID(),
                         fwdtrack.chi2MatchMCHMFT(),
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

  template <typename TBC>
  void fillFIThistograms(TBC const& bc)
  {
    LOGF(debug, "");
    std::array<bool, 5> triggers{{true, !udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits()),
                                  udhelpers::TVX(bc), udhelpers::TSC(bc), udhelpers::TCE(bc)}};
    LOGF(debug, "triggers %d %d %d %d %d", triggers[0], triggers[1], triggers[2], triggers[3], triggers[4]);
    if (bc.has_foundFV0()) {
      auto fv0 = bc.foundFV0();
      auto ampA = udhelpers::FV0AmplitudeA(fv0);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, 0);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, triggers[1] ? 1 : 5);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, triggers[2] ? 2 : 6);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, triggers[3] ? 3 : 7);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, triggers[4] ? 4 : 8);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kIsBBV0A) ? 9 : 10);
      registry.get<TH2>(HIST("reco/fv0"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kNoBGV0A) ? 11 : 12);
    }
    if (bc.has_foundFT0()) {
      auto ft0 = bc.foundFT0();
      auto ampA = udhelpers::FT0AmplitudeA(ft0);
      auto ampC = udhelpers::FT0AmplitudeC(ft0);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, 0);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, 0);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, triggers[1] ? 1 : 5);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, triggers[1] ? 1 : 5);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, triggers[2] ? 2 : 6);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, triggers[2] ? 2 : 6);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, triggers[3] ? 3 : 7);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, triggers[3] ? 3 : 7);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, triggers[4] ? 4 : 8);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, triggers[4] ? 4 : 8);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kIsBBT0A) ? 9 : 10);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, bc.selection_bit(o2::aod::evsel::kIsBBT0C) ? 9 : 10);
      registry.get<TH2>(HIST("reco/ft0A"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kNoBGT0A) ? 11 : 12);
      registry.get<TH2>(HIST("reco/ft0C"))->Fill(ampC, bc.selection_bit(o2::aod::evsel::kNoBGT0C) ? 11 : 12);
    }
    if (bc.has_foundFDD()) {
      auto fdd = bc.foundFDD();
      auto ampA = udhelpers::FDDAmplitudeA(fdd);
      auto ampC = udhelpers::FDDAmplitudeC(fdd);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, 0);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, 0);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, triggers[1] ? 1 : 5);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, triggers[1] ? 1 : 5);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, triggers[2] ? 2 : 6);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, triggers[2] ? 2 : 6);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, triggers[3] ? 3 : 7);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, triggers[3] ? 3 : 7);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, triggers[4] ? 4 : 8);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, triggers[4] ? 4 : 8);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kIsBBFDA) ? 9 : 10);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, bc.selection_bit(o2::aod::evsel::kIsBBFDC) ? 9 : 10);
      registry.get<TH2>(HIST("reco/fddA"))->Fill(ampA, bc.selection_bit(o2::aod::evsel::kNoBGFDA) ? 11 : 12);
      registry.get<TH2>(HIST("reco/fddC"))->Fill(ampC, bc.selection_bit(o2::aod::evsel::kNoBGFDC) ? 11 : 12);
    }
  }

  void init(InitContext&)
  {
    LOGF(debug, "<DGCandProducer> beginning of init reached");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
    diffCuts = (DGCutparHolder)DGCuts;

    const int nXbinsInStatH = 25;

    // add histograms for the different process functions
    registry.add("reco/Stat", "Cut statistics;; Collisions", {HistType::kTH1F, {{nXbinsInStatH, -0.5, static_cast<float>(nXbinsInStatH - 0.5)}}});
    registry.add("reco/pt1Vspt2", "2 prong events, p_{T} versus p_{T}", {HistType::kTH2F, {{100, -3., 3.}, {100, -3., 3.0}}});
    registry.add("reco/TPCsignal1", "2 prong events, TPC signal versus p_{T} of particle 1", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
    registry.add("reco/TPCsignal2", "2 prong events, TPC signal versus p_{T} of particle 2", {HistType::kTH2F, {{200, -3., 3.}, {200, 0., 100.0}}});
    registry.add("reco/sig1VsSig2TPC", "2 prong events, TPC signal versus TPC signal", {HistType::kTH2F, {{100, 0., 100.}, {100, 0., 100.}}});

    // FIT amplitudes
    //   0: unconditional
    //   1: TOR              5: no TOR
    //   2: TVX              6: no TVX
    //   3: TSC              7: no TSC
    //   4: TCE              8: no TCE
    //   9: IsBBXXX         10: !IsBBXXX
    //  11: kNoBGXXX        12: !kNoBGXXX
    const int nXbinsFITH = 201;
    registry.add("reco/fv0", "FV0 amplitudes", {HistType::kTH2F, {{nXbinsFITH, -0.5, nXbinsFITH - 0.5}, {13, -0.5, 12.5}}});
    registry.add("reco/ft0A", "FT0A amplitudes", {HistType::kTH2F, {{nXbinsFITH, -0.5, nXbinsFITH - 0.5}, {13, -0.5, 12.5}}});
    registry.add("reco/ft0C", "FT0C amplitudes", {HistType::kTH2F, {{nXbinsFITH, -0.5, nXbinsFITH - 0.5}, {13, -0.5, 12.5}}});
    registry.add("reco/fddA", "FDDA amplitudes", {HistType::kTH2F, {{nXbinsFITH, -0.5, nXbinsFITH - 0.5}, {13, -0.5, 12.5}}});
    registry.add("reco/fddC", "FDDC amplitudes", {HistType::kTH2F, {{nXbinsFITH, -0.5, nXbinsFITH - 0.5}, {13, -0.5, 12.5}}});

    std::string labels[nXbinsInStatH] = {"all", "hasBC", "accepted", "FITveto", "MID trk", "global not PV trk", "not global PV trk",
                                         "ITS-only PV trk", "TOF PV trk fraction", "n PV trks", "PID", "pt", "eta", "net charge",
                                         "inv mass", "evsel TF border", "evsel no pile-up", "evsel ITSROF", "evsel z-vtx", "evsel ITSTPC vtx", "evsel TRD vtx", "evsel TOF vtx", "", "", ""};

    registry.get<TH1>(HIST("reco/Stat"))->SetNdivisions(nXbinsInStatH, "X");
    for (int iXbin(1); iXbin < nXbinsInStatH + 1; iXbin++) {
      registry.get<TH1>(HIST("reco/Stat"))->GetXaxis()->ChangeLabel(iXbin, 45, 0.03, 33, -1, -1, labels[iXbin - 1]);
    }

    LOGF(debug, "<DGCandProducer> end of init reached");
  }

  // process function for real data
  void process(CC const& collision, BCs const& bcs, TCs& tracks, FWs& fwdtracks,
               aod::Zdcs& /*zdcs*/, aod::FV0As& fv0as, aod::FT0s& ft0s, aod::FDDs& fdds)
  {
    LOGF(debug, "<DGCandProducer>  collision %d", collision.globalIndex());
    registry.get<TH1>(HIST("reco/Stat"))->Fill(0., 1.);

    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    auto bc = collision.foundBC_as<BCs>();
    int trs = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      trs = 1;
    }
    int trofs = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      trofs = 1;
    }
    int hmpr = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      hmpr = 1;
    }
    double ir = 0.;
    const uint64_t ts = bc.timestamp();
    const int runnumber = bc.runNumber();
    if (bc.has_zdc()) {
      ir = mRateFetcher.fetch(ccdb.service, ts, runnumber, "ZNC hadronic") * 1.e-3;
    }
    uint8_t chFT0A = 0;
    uint8_t chFT0C = 0;
    uint8_t chFDDA = 0;
    uint8_t chFDDC = 0;
    uint8_t chFV0A = 0;
    int occ = 0;
    occ = collision.trackOccupancyInTimeRange();
    LOGF(debug, "<DGCandProducer>  BC id %d", bc.globalBC());

    // fill FIT histograms
    fillFIThistograms(bc);

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
    LOGF(debug, "<DGCandProducer>  Size of bcRange %d", bcRange.size());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, tracks, fwdtracks);

    // save DG candidates
    registry.get<TH1>(HIST("reco/Stat"))->Fill(isDGEvent + 2, 1.);
    if (isDGEvent == 0) {
      LOGF(debug, "<DGCandProducer>  Data: good collision!");

      // fill FITInfo
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bc, bcs, ft0s, fv0as, fdds);

      // update DG candidates tables
      auto rtrwTOF = udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib());
      int upc_flag = 0;
      ushort flags = collision.flags();
      if (flags & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode)
        upc_flag = 1;
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(), upc_flag,
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
      outputCollisionSelExtras(chFT0A, chFT0C, chFDDA, chFDDC, chFV0A, occ, ir, trs, trofs, hmpr);
      outputCollsLabels(collision.globalIndex());

      // update DGTracks tables
      for (auto& track : tracks) {
        if (saveAllTracks || track.isPVContributor()) {
          updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
        }
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

        float timeZNA = zdc.timeZNA();
        float timeZNC = zdc.timeZNC();
        float eComZNA = zdc.energyCommonZNA();
        float eComZNC = zdc.energyCommonZNC();
        outputZdcsReduced(outputCollisions.lastIndex(), timeZNA, timeZNC, eComZNA, eComZNC);
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
    auto bc = mccol.template bc_as<BCs>();
    outputMcCollisions(bc.globalBC(),
                       mccol.generatorsID(),
                       mccol.posX(),
                       mccol.posY(),
                       mccol.posZ(),
                       mccol.t(),
                       mccol.weight(),
                       mccol.impactParameter());
  }

  template <typename TMcParticle>
  void updateUDMcParticle(TMcParticle const& McPart, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // save McPart
    // mother and daughter indices are set to -1
    // ATTENTION: this can be improved to also include mother and daughter indices
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};

    // update UDMcParticles
    if (mcPartIsSaved.find(McPart.globalIndex()) == mcPartIsSaved.end()) {
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
      mcPartIsSaved[McPart.globalIndex()] = outputMcParticles.lastIndex();
    }
  }

  template <typename TMcParticles>
  void updateUDMcParticles(TMcParticles const& McParts, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    LOGF(debug, "number of McParticles %d", McParts.size());
    // save McParts
    // new mother and daughter ids
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};
    int64_t newval = -1;

    // Determine the particle indices within the UDMcParticles table
    // before filling the table
    // This is needed to be able to assign the new daughter indices
    std::map<int64_t, int64_t> oldnew;
    auto lastId = outputMcParticles.lastIndex();
    for (auto mcpart : McParts) {
      auto oldId = mcpart.globalIndex();
      if (mcPartIsSaved.find(oldId) != mcPartIsSaved.end()) {
        oldnew[oldId] = mcPartIsSaved[oldId];
      } else {
        lastId++;
        oldnew[oldId] = lastId;
      }
    }

    // all particles of the McCollision are saved
    for (auto mcpart : McParts) {
      LOGF(debug, "  p (%d) %d", mcpart.pdgCode(), mcpart.globalIndex());
      if (mcPartIsSaved.find(mcpart.globalIndex()) == mcPartIsSaved.end()) {
        // mothers
        newmids.clear();
        auto oldmids = mcpart.mothersIds();
        for (auto oldmid : oldmids) {
          auto m = McParts.rawIteratorAt(oldmid);
          LOGF(debug, "    m %d", m.globalIndex());
          if (mcPartIsSaved.find(oldmid) != mcPartIsSaved.end()) {
            newval = mcPartIsSaved[oldmid];
          } else {
            newval = -1;
          }
          LOGF(debug, "    mid o %i n %i", oldmid, newval);
          newmids.push_back(newval);
        }

        // daughters
        auto olddids = mcpart.daughtersIds();
        for (uint ii = 0; ii < olddids.size(); ii++) {
          if (oldnew.find(olddids[ii]) != oldnew.end()) {
            newval = oldnew[olddids[ii]];
          } else {
            newval = -1;
          }
          LOGF(debug, "    did o %i n %i", olddids[ii], newval);
          newdids[ii] = newval;
        }
        LOGF(debug, "    ms %i ds %i", oldmids.size(), olddids.size());

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
        mcPartIsSaved[mcpart.globalIndex()] = outputMcParticles.lastIndex();
        LOGF(debug, "  mcpart %d -> udmcpart %d", mcpart.globalIndex(), mcPartIsSaved[mcpart.globalIndex()]);
      }
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabel(TTrack const& udtrack, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
    auto trackId = udtrack.trackId();
    if (trackId >= 0) {
      auto track = udtrack.template track_as<TCs>();
      auto mcTrackId = track.mcParticleId();
      if (mcTrackId >= 0) {
        if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
          outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, track.mcMask());
      }
    } else {
      outputMcTrackLabels(-1, -1);
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabels(TTrack const& udtracks, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // loop over all tracks
    for (auto udtrack : udtracks) {
      // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
      auto trackId = udtrack.trackId();
      if (trackId >= 0) {
        auto track = udtrack.template track_as<TCs>();
        auto mcTrackId = track.mcParticleId();
        if (mcTrackId >= 0) {
          if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
            outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
          } else {
            outputMcTrackLabels(-1, track.mcMask());
          }
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
                 CCs const& /*collisions*/, BCs const& /*bcs*/, TCs const& /*tracks*/)
  {
    LOGF(info, "Number of McCollisions %d", mccols.size());
    LOGF(info, "Number of DG candidates %d", dgcands.size());
    LOGF(info, "Number of UD tracks %d", udtracks.size());
    if (dgcands.size() <= 0) {
      LOGF(info, "No DG candidates to save!");
      return;
    }

    // use a hash table to keep track of the McCollisions which have been added to the UDMcCollision table
    // {McCollisionId : udMcCollisionId}
    // similar for the McParticles which have been added to the UDMcParticle table
    // {McParticleId : udMcParticleId}
    std::map<int64_t, int64_t> mcColIsSaved;
    std::map<int64_t, int64_t> mcPartIsSaved;

    // loop over McCollisions and UDCCs simultaneously
    auto mccol = mccols.iteratorAt(0);
    auto dgcand = dgcands.iteratorAt(0);
    auto lastmccol = mccols.iteratorAt(mccols.size() - 1);
    auto lastdgcand = dgcands.iteratorAt(dgcands.size() - 1);

    // advance dgcand and mccol until both are AtEnd
    int64_t mccolId = mccol.globalIndex();
    int64_t mcdgId = -1;
    int64_t colId = -1;
    auto dgcandAtEnd = dgcand == lastdgcand;
    auto mccolAtEnd = mccol == lastmccol;
    bool goon = !dgcandAtEnd || !mccolAtEnd;

    while (goon) {
      // check if dgcand has an associated Collision and McCollision
      if (dgcand.has_collision()) {
        auto dgcandCol = dgcand.collision_as<CCs>();
        colId = dgcandCol.globalIndex();
        if (dgcandCol.has_mcCollision()) {
          mcdgId = dgcandCol.mcCollision().globalIndex();
        } else {
          mcdgId = -1;
        }
      } else {
        colId = -1;
        mcdgId = -1;
      }
      LOGF(debug, "");
      LOGF(debug, "dgcand %d mcdgId %d colId %d mccolId %d - UDMcCollsLabels %d UDMcCollisions %d", dgcand.globalIndex(), mcdgId, colId, mccolId, outputMcCollsLabels.lastIndex(), outputMcCollisions.lastIndex());

      // two cases to consider
      // 1. mcdgId <= mccolId: the event to process is a dgcand. In this case the Mc tables as well as the McLabel tables are updated
      // 2. mccolId < mcdgId: the event to process is an MC event of interest without reconstructed dgcand. In this case only the Mc tables are updated
      if ((!dgcandAtEnd && !mccolAtEnd && (mcdgId <= mccolId)) || mccolAtEnd) {
        // this is case 1.

        // update UDMcCollisions and UDMcColsLabels (for each UDCollision -> UDMcCollisions)
        // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
        // get dgcand tracks
        auto dgTracks = udtracks.sliceByCached(aod::udtrack::udCollisionId, dgcand.globalIndex(), cache);

        // If the dgcand has an associated McCollision then the McCollision and all associated
        // McParticles are saved
        if (mcdgId >= 0) {
          if (mcColIsSaved.find(mcdgId) == mcColIsSaved.end()) {
            // update UDMcCollisions
            LOGF(debug, "  writing mcCollision %d to UDMcCollisions", mcdgId);
            auto dgcandMcCol = dgcand.collision_as<CCs>().mcCollision();
            updateUDMcCollisions(dgcandMcCol);
            mcColIsSaved[mcdgId] = outputMcCollisions.lastIndex();
          }

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          LOGF(debug, "  writing %d to outputMcCollsLabels", mcColIsSaved[mcdgId]);
          outputMcCollsLabels(mcColIsSaved[mcdgId]);

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mcdgId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mcdgId], mcPartIsSaved);

          // update UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          updateUDMcTrackLabels(dgTracks, mcPartIsSaved);

        } else {
          // If the dgcand has no associated McCollision then only the McParticles which are associated
          // with the tracks of the dgcand are saved

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          LOGF(debug, "  writing %d to UDMcCollsLabels", -1);
          outputMcCollsLabels(-1);

          // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          // loop over tracks of dgcand
          for (auto dgtrack : dgTracks) {
            if (dgtrack.has_track()) {
              auto track = dgtrack.track_as<TCs>();
              if (track.has_mcParticle()) {
                auto mcPart = track.mcParticle();
                updateUDMcParticle(mcPart, -1, mcPartIsSaved);
                updateUDMcTrackLabel(dgtrack, mcPartIsSaved);
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

        // update UDMcCollisions and UDMcParticles
        if (mcColIsSaved.find(mccolId) == mcColIsSaved.end()) {

          // update UDMcCollisions
          LOGF(debug, "  writing mcCollision %d to UDMcCollisions", mccolId);
          updateUDMcCollisions(mccol);
          mcColIsSaved[mccolId] = outputMcCollisions.lastIndex();

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccolId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mccolId], mcPartIsSaved);
        }

        // advance mccol
        if (mccol != lastmccol) {
          mccol++;
          mccolId = mccol.globalIndex();
        } else {
          mccolAtEnd = true;
        }
      }
      LOGF(debug, "  UDMcCollsLabels %d (of %d) UDMcCollisions %d", outputMcCollsLabels.lastIndex(), dgcands.size() - 1, outputMcCollisions.lastIndex());
      goon = !dgcandAtEnd || !mccolAtEnd;
    }
  }
  PROCESS_SWITCH(McDGCandProducer, processMC, "Produce MC tables", false);

  void processDummy(aod::Collisions const& /*collisions*/)
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
