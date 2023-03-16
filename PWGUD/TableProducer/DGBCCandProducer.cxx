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
/// \brief
///
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  30.09.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/DGSelector.h"
#include "DGBCCandProducer.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// -----------------------------------------------------------------------------
struct tracksWGTInBCs {
  Produces<aod::TracksWGTInBCs> tracksWGTInBCs;
  Produces<aod::FwdTracksWGTInBCs> fwdTracksWGTInBCs;

  HistogramRegistry registry{
    "registry",
    {}};

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using TC = TCs::iterator;
  using ATs = aod::AmbiguousTracks;

  Preslice<aod::AmbiguousTracks> perTrack = aod::ambiguous::trackId;
  Preslice<aod::AmbiguousFwdTracks> perFwdTrack = aod::ambiguous::fwdtrackId;

  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processBarrel")) {
      registry.add("barrelTracks", "#barrelTracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    }
    if (context.mOptions.get<bool>("processForward")) {
      registry.add("forwardTracks", "#forwardTracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    }
  }

  void processBarrel(BCs const& bcs, CCs const& collisions, TCs const& tracks, ATs const& ambTracks)
  {
    // run number
    int rnum = bcs.iteratorAt(0).runNumber();

    // container to sort tracks with good timing according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> tracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all tracks and fill tracksInBCList
    LOGF(info, "Number of barrel tracks: %d", tracks.size());
    for (auto const& track : tracks) {
      registry.get<TH1>(HIST("barrelTracks"))->Fill(0., 1.);
      auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
      if (ambTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(2., 1.);
      } else {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(1., 1.);
      }

      // only consider tracks with good timing
      if (track.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(3., 1.);

        // get first compatible BC
        // auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
        if (ambTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("barrelTracks"))->Fill(4., 1.);

          // compute the BC closest in time
          auto firstCompatibleBC = ambTracksSlice.begin().bc().begin().globalBC();
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        } else {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          closestBC = track.collision_as<CCs>().foundBC_as<BCs>().globalBC();
        }

        // update tracksInBCList
        tracksInBCList[closestBC].emplace_back((int32_t)track.globalIndex());
      }
    }

    // fill tracksWGTInBCs
    int indBCToStart = 0;
    int indBCToSave;
    for (auto const& tracksInBC : tracksInBCList) {
      indBCToSave = -1;
      if (tracksInBC.second.size() > 0) {
        // find corresponding BC
        for (auto ind = indBCToStart; ind < bcs.size(); ind++) {
          auto bc = bcs.rawIteratorAt(ind);
          if (bc.globalBC() == tracksInBC.first) {
            indBCToSave = ind;
            indBCToStart = ind;
            break;
          }
          if (bc.globalBC() > tracksInBC.first) {
            break;
          }
        }
        tracksWGTInBCs(indBCToSave, rnum, tracksInBC.first, tracksInBC.second);
        LOGF(debug, " BC %i/%u with %i tracks with good timing", indBCToSave, tracksInBC.first, tracksInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWGTInBCs, processBarrel, "Process barrel tracks", false);

  void processForward(BCs& bcs, CCs& collisions, aod::FwdTracks& fwdTracks, aod::AmbiguousFwdTracks& ambFwdTracks)
  {
    // run number
    int rnum = bcs.iteratorAt(0).runNumber();

    // container to sort forward tracks according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> fwdTracksWGTInBCList{};
    uint64_t closestBC = 0;

    // loop over all forward tracks and fill fwdTracksWGTInBCList
    LOGF(info, "Number of forward tracks: %d", fwdTracks.size());
    for (auto const& fwdTrack : fwdTracks) {
      registry.get<TH1>(HIST("forwardTracks"))->Fill(0., 1.);
      auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
      if (ambFwdTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(2., 1.);
      } else {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(1., 1.);
      }

      // only consider tracks with trackTimeRes < LHCBunchSpacingNS
      LOGF(debug, "Time resolution of fwdTrack %f", fwdTrack.trackTimeRes());
      if (fwdTrack.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(3., 1.);

        // get first compatible BC
        // auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
        if (ambFwdTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("forwardTracks"))->Fill(4., 1.);

          // compute the BC closest in time
          auto firstCompatibleBC = ambFwdTracksSlice.begin().bc().begin().globalBC();
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        } else {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          closestBC = fwdTrack.collision_as<CCs>().bc_as<BCs>().globalBC();
        }

        // update tracksInBCList
        fwdTracksWGTInBCList[closestBC].emplace_back((int32_t)fwdTrack.globalIndex());
      }
    }

    // fill fwdTracksWGTInBCs
    int indBCToStart = 0;
    int indBCToSave;
    for (auto const& fwdTracksWGTInBC : fwdTracksWGTInBCList) {
      indBCToSave = -1;
      if (fwdTracksWGTInBC.second.size() > 0) {
        // find corresponding BC
        for (auto ind = indBCToStart; ind < bcs.size(); ind++) {
          auto bc = bcs.rawIteratorAt(ind);
          if (bc.globalBC() == fwdTracksWGTInBC.first) {
            indBCToSave = ind;
            indBCToStart = ind;
            break;
          }
          if (bc.globalBC() > fwdTracksWGTInBC.first) {
            break;
          }
        }
        fwdTracksWGTInBCs(indBCToSave, rnum, fwdTracksWGTInBC.first, fwdTracksWGTInBC.second);
        LOGF(debug, " BC %i/%u with %i forward tracks with good timing", indBCToSave, fwdTracksWGTInBC.first, fwdTracksWGTInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWGTInBCs, processForward, "Process forward tracks", false);

  void processNone(BCs&)
  {
    LOGF(info, "Tables tibcs and ftibcs are not filled!");
  }
  PROCESS_SWITCH(tracksWGTInBCs, processNone, "Process dummy task", true);
};

// -----------------------------------------------------------------------------
struct DGBCCandProducer {
  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDTracks> outputTracks;
  // Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector = DGSelector();

  HistogramRegistry registry{
    "registry",
    {}};

  using TIBCs = aod::TracksWGTInBCs;
  using TIBC = TIBCs::iterator;
  using FTIBCs = aod::FwdTracksWGTInBCs;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FTCs = aod::FwdTracks;

  Preslice<CCs> CCperBC = aod::evsel::foundBCId;
  Preslice<TCs> TCperCollision = aod::track::collisionId;
  Preslice<aod::FwdTracks> FWperCollision = aod::fwdtrack::collisionId;
  Preslice<TIBCs> TIBCperBC = aod::dgbcandidate::bcId;
  Preslice<FTIBCs> FTIBCperBC = aod::dgbcandidate::bcId;

  // fill BB and BG information into FITInfo
  template <typename BCR>
  void fillBGBBFlags(upchelpers::FITInfo& info, uint64_t const& minbc, BCR const& bcrange)
  {
    for (auto const& bc2u : bcrange) {

      // 0 <= bit <= 31
      auto bit = bc2u.globalBC() - minbc;
      if (!bc2u.selection()[evsel::kNoBGT0A])
        SETBIT(info.BGFT0Apf, bit);
      if (!bc2u.selection()[evsel::kNoBGT0C])
        SETBIT(info.BGFT0Cpf, bit);
      if (bc2u.selection()[evsel::kIsBBT0A])
        SETBIT(info.BBFT0Apf, bit);
      if (bc2u.selection()[evsel::kIsBBT0C])
        SETBIT(info.BBFT0Cpf, bit);
      if (!bc2u.selection()[evsel::kNoBGV0A])
        SETBIT(info.BGFV0Apf, bit);
      if (bc2u.selection()[evsel::kIsBBV0A])
        SETBIT(info.BBFV0Apf, bit);
      if (!bc2u.selection()[evsel::kNoBGFDA])
        SETBIT(info.BGFDDApf, bit);
      if (!bc2u.selection()[evsel::kNoBGFDC])
        SETBIT(info.BGFDDCpf, bit);
      if (bc2u.selection()[evsel::kIsBBFDA])
        SETBIT(info.BBFDDApf, bit);
      if (bc2u.selection()[evsel::kIsBBFDC])
        SETBIT(info.BBFDDCpf, bit);
    }
  }

  // extract FIT information
  upchelpers::FITInfo getFITinfo(uint64_t const& bcnum, BCs const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // FITinfo
    upchelpers::FITInfo info{};
    uint64_t minbc = bcnum > 16 ? bcnum - 16 : 0;

    // find bc with globalBC = bcnum
    Partition<BCs> selbc = aod::bc::globalBC == bcnum;
    selbc.bindTable(bcs);

    // if BC exists then update FIT information for this BC
    if (selbc.size() > 0) {
      auto bc = bcs.iteratorAt(selbc.begin().globalIndex());

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

      auto bcrange = udhelpers::compatibleBCs(bc, bcnum, 16, bcs);
      fillBGBBFlags(info, minbc, bcrange);
    } else {
      auto bcrange = udhelpers::compatibleBCs(bcnum, 16, bcs);
      fillBGBBFlags(info, minbc, bcrange);
    }
    return info;
  }

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack>
  void updateUDTrackTables(TTrack const& track, uint64_t const& bcnum)
  {
    outputTracks(outputCollisions.lastIndex(),
                 track.px(), track.py(), track.pz(), track.sign(),
                 bcnum, track.trackTime(), track.trackTimeRes());
    // outputTracksCov(track.x(), track.y(), track.z(), track.sigmaY(), track.sigmaZ());
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

  // update UDTables
  template <typename TTracks>
  void updateUDTables(bool onlyPV, uint64_t bcnum, int rnum, float vx, float vy, float vz,
                      uint16_t const& ntrks, int8_t const& ncharge, float const& rtrwTOF,
                      TTracks const& tracks, upchelpers::FITInfo const& fitInfo)
  {
    outputCollisions(bcnum, rnum, vx, vy, vz, ntrks, ncharge, rtrwTOF);
    outputCollisionsSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C,
                         fitInfo.triggerMaskFT0,
                         fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC,
                         fitInfo.triggerMaskFDD,
                         fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                         fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                         fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                         fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);

    // update DGTracks tables
    for (auto const& track : tracks) {
      if (track.isPVContributor() || !onlyPV) {
        updateUDTrackTables(track, bcnum);
      }
    }
  }

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    if (context.mOptions.get<bool>("processTable") || context.mOptions.get<bool>("processData")) {
      registry.add("bcFlag", "#bcFlag", {HistType::kTH1F, {{64, -0.5, 63.5}}});
      registry.add("candCase", "#candCase", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.add("isDG1vsisDG2", "#isDG1vsisDG2", {HistType::kTH2F, {{13, -1.5, 11.5}, {13, -1.5, 11.5}}});
      registry.add("ntr1vsntr2All", "#ntr1vsntr2All", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("ntr1vsntr2Cand", "#ntr1vsntr2Cand", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("ptvsdcaxy1", "#ptvsdcaxy1", {HistType::kTH2F, {{50, 0., 5.}, {50, -5., 5.}}});
      registry.add("ptvsdcaxy2", "#ptvsdcaxy2", {HistType::kTH2F, {{50, 0., 5.}, {50, -5., 5.}}});
    }
  }

  void processTable(TIBC const& tibc, BCs const& bcs, CCs const& collisions,
                    TCs const& tracks, aod::FwdTracks const& fwdtracks, FTIBCs const& ftibcs,
                    aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // fill FITInfo
    auto bcnum = tibc.bcnum();
    upchelpers::FITInfo fitInfo = getFITinfo(bcnum, bcs, ft0s, fv0as, fdds);

    // check if DG event
    // distinguish between cases with and without associated BC
    // 1. candidate has associated BC and associated collision
    // 2. candidate has associated BC but no associated collision
    // 3. candidate has no associated BC
    int isDG = -1;
    float rtrwTOF = -1.;
    int8_t nCharge;
    if (tibc.has_bc()) {
      LOGF(debug, "[1.,2.] BC found");

      // get associated bc
      auto bc = tibc.bc_as<BCs>();

      // is there an associated collision?
      auto colSlize = collisions.sliceBy(CCperBC, bc.globalIndex());
      if (colSlize.size() > 0) {
        LOGF(debug, "  1. BC has collision");
        auto col = colSlize.rawIteratorAt(0);
        auto colTracks = tracks.sliceBy(TCperCollision, col.globalIndex());
        auto colFwdTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
        auto bcRange = udhelpers::compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
        isDG = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);

        // update UDTables
        if (isDG == 0) {
          registry.get<TH1>(HIST("candCase"))->Fill(1, 1.);
          rtrwTOF = udhelpers::rPVtrwTOF<true>(colTracks, col.numContrib());
          nCharge = udhelpers::netCharge<true>(colTracks);

          outputCollisions(bc.globalBC(), bc.runNumber(),
                           col.posX(), col.posY(), col.posZ(),
                           col.numContrib(), nCharge,
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
          for (auto const& track : colTracks) {
            updateUDTrackTables(track, bc.globalBC());
          }
        }
      } else {
        LOGF(debug, "  2. BC has NO collision");
        auto tracksArray = tibc.track_as<TCs>();
        auto bcRange = udhelpers::compatibleBCs(bc, bc.globalBC(), diffCuts.minNBCs(), bcs);

        // does BC have fwdTracks?
        if (ftibcs.size() > 0) {
          auto ftibcSlice = ftibcs.sliceBy(FTIBCperBC, bc.globalIndex());
          if (ftibcSlice.size() > 0) {
            auto fwdTracksArray = ftibcSlice.begin().fwdtrack_as<FTCs>();
            isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
          } else {
            auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
            isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
          }
        } else {
          auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
          isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        }

        // update UDTables
        if (isDG == 0) {
          registry.get<TH1>(HIST("candCase"))->Fill(2, 1.);
          rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
          nCharge = udhelpers::netCharge<false>(tracksArray);

          outputCollisions(bc.globalBC(), bc.runNumber(),
                           -1., 1., -1.,
                           tracksArray.size(), nCharge,
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
          for (auto const& track : tracksArray) {
            updateUDTrackTables(track, bc.globalBC());
          }
        }
      }
    } else {
      LOGF(debug, "3. BC NOT found");

      // the BC is not contained in the BCs table
      auto tracksArray = tibc.track_as<TCs>();
      auto bcRange = udhelpers::compatibleBCs(bcnum, diffCuts.minNBCs(), bcs);

      // does BC have fwdTracks?
      if (ftibcs.size() > 0) {
        Partition<FTIBCs> ftibcPart = aod::dgbcandidate::bcnum == bcnum;
        ftibcPart.bindTable(ftibcs);

        if (ftibcPart.size() > 0) {
          auto fwdTracksArray = ftibcPart.begin().fwdtrack_as<FTCs>();
          isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        } else {
          auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
          isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        }
      } else {
        auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
        isDG = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
      }

      // update UDTables
      if (isDG == 0) {
        registry.get<TH1>(HIST("candCase"))->Fill(3, 1.);
        rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
        nCharge = udhelpers::netCharge<false>(tracksArray);

        outputCollisions(bcnum, tibc.runNumber(),
                         -2., 2., -2.,
                         tracksArray.size(), nCharge,
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
        for (auto const& track : tracksArray) {
          updateUDTrackTables(track, bcnum);
        }
      }
    }
  }

  PROCESS_SWITCH(DGBCCandProducer, processTable, "Produce UDTables", true);

  void processData(BCs const& bcs, CCs const& collisions,
                   TCs const& tracks, FTCs const& fwdtracks, TIBCs const& tibcs, FTIBCs const& ftibcs,
                   aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    int isDG1, isDG2;
    int ntr1, ntr2;
    // flag BCs
    //  0: all
    //  1: in BCs table
    //  2: has associated collision
    //  3: is DG1 candidate
    //  4: in TIBCs table
    //  5: is DG2 candidate
    int8_t bcFlag = 0;

    // loop over all possible BCs
    if (bcs.size() <= 0) {
      return;
    }

    // run over globalBC [minGlobalBC, maxGlobalBC] ...
    uint64_t minGlobalBC = bcs.iteratorAt(0).globalBC();
    uint64_t maxGlobalBC = bcs.iteratorAt(bcs.size() - 1).globalBC();
    // ... and advance these pointers step-by-step
    auto bc = bcs.iteratorAt(0);
    auto tibc = tibcs.iteratorAt(0);
    auto ftibc = ftibcs.iteratorAt(0);
    auto lastbc = bcs.iteratorAt(bcs.size() - 1);
    auto lasttibc = tibcs.iteratorAt(tibcs.size() - 1);
    auto lastftibc = ftibcs.iteratorAt(ftibcs.size() - 1);
    LOGF(info, "bcs %d tibcs %d ftibcs %d", bcs.size(), tibcs.size(), ftibcs.size());

    for (auto bcnum = minGlobalBC; bcnum <= maxGlobalBC; bcnum++) {
      LOGF(debug, "max %d now %d", maxGlobalBC, bcnum);

      // reset counters
      bcFlag = 1; // bit 0 is always set
      isDG1 = -1;
      isDG2 = -1;
      ntr1 = -1;
      ntr2 = -1;

      // find BC in BCs table
      while (bc.globalBC() < bcnum && bc != lastbc) {
        ++bc;
      }
      if (bc.globalBC() == bcnum) {
        SETBIT(bcFlag, 1);

        // find associated collision
        auto colSlize = collisions.sliceBy(CCperBC, bc.globalIndex());
        if (colSlize.size() > 0) {
          SETBIT(bcFlag, 2);

          auto col = colSlize.begin();
          ntr1 = col.numContrib();
          auto colTracks = tracks.sliceBy(TCperCollision, col.globalIndex());
          auto colFwdTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
          auto bcRange = udhelpers::compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
          isDG1 = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);
          if (isDG1 == 0) {
            // this is a DG candidate with proper collision vertex
            SETBIT(bcFlag, 3);
            registry.get<TH1>(HIST("candCase"))->Fill(1, 1.);

            auto rtrwTOF = udhelpers::rPVtrwTOF<true>(colTracks, col.numContrib());
            auto nCharge = udhelpers::netCharge<true>(colTracks);
            auto fitInfo = getFITinfo(bcnum, bcs, ft0s, fv0as, fdds);
            updateUDTables(false, bcnum, bc.runNumber(), col.posX(), col.posY(), col.posZ(),
                           col.numContrib(), nCharge, rtrwTOF, colTracks, fitInfo);
          }
        }
      }

      // find BC in TIBCs table
      if (ftibcs.size() > 0) {
        while (tibc.bcnum() < bcnum && tibc != lasttibc) {
          ++tibc;
        }
        if (tibc.bcnum() == bcnum) {
          SETBIT(bcFlag, 4);

          auto bcRange = udhelpers::compatibleBCs(bc, bcnum, diffCuts.minNBCs(), bcs);
          auto tracksArray = tibc.track_as<TCs>();
          ntr2 = tracksArray.size();

          // find BC in FTIBCs table
          if (ftibcs.size() > 0) {
            while (ftibc.bcnum() < bcnum && ftibc != lastftibc) {
              ++ftibc;
            }
            if (ftibc.bcnum() == bcnum) {
              auto fwdTracksArray = ftibc.fwdtrack_as<FTCs>();
              isDG2 = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
            } else {
              auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
              isDG2 = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
            }
          } else {
            auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
            isDG2 = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
          }

          if (isDG2 == 0) {
            // this is a DG candidate with tracks-in-BC
            SETBIT(bcFlag, 5);
            registry.get<TH1>(HIST("candCase"))->Fill(2, 1.);

            auto rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
            auto nCharge = udhelpers::netCharge<false>(tracksArray);
            auto fitInfo = getFITinfo(bcnum, bcs, ft0s, fv0as, fdds);
            updateUDTables(false, bcnum, tibc.runNumber(), -1., 1., -1,
                           tracksArray.size(), nCharge, rtrwTOF, tracksArray, fitInfo);
          }
        }
      }

      // update histograms
      registry.get<TH1>(HIST("bcFlag"))->Fill(bcFlag, 1.);
      registry.get<TH2>(HIST("isDG1vsisDG2"))->Fill(isDG1, isDG2);
      registry.get<TH2>(HIST("ntr1vsntr2All"))->Fill(ntr1, ntr2);
      if (isDG1 == 0 && isDG2 == 0) {
        registry.get<TH2>(HIST("ntr1vsntr2Cand"))->Fill(ntr1, ntr2);
      }
    }
  }

  PROCESS_SWITCH(DGBCCandProducer, processData, "Produce UDTables", true);

  void processDummy(BCs const& bcs, CCs const& collisions,
                    TCs const& tracks, FTCs const& fwdtracks, TIBCs const& tibcs, FTIBCs const& ftibcs,
                    aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    LOGF(info, "size of");
    LOGF(info, "bcs %d", bcs.size());
    LOGF(info, "collisions %d", collisions.size());
    LOGF(info, "tracks %d", tracks.size());
    LOGF(info, "fwdtracks %d", fwdtracks.size());
  }
  PROCESS_SWITCH(DGBCCandProducer, processDummy, "Dummy task", true);
};

// -----------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tracksWGTInBCs>(cfgc, TaskName{"trackswgtinbcs"}),
    adaptAnalysisTask<DGBCCandProducer>(cfgc, TaskName{"dgbccandproducer"}),
  };
}

// -----------------------------------------------------------------------------
