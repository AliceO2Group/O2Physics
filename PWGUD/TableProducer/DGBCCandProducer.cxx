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

#include <algorithm>
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
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  using TC = TCs::iterator;
  using ATs = aod::AmbiguousTracks;

  Preslice<aod::AmbiguousTracks> perTrack = aod::ambiguous::trackId;
  Preslice<aod::AmbiguousFwdTracks> perFwdTrack = aod::ambiguous::fwdtrackId;

  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processBarrel")) {
      registry.add("barrel/Tracks", "Number of barrel track types", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    }
    if (context.mOptions.get<bool>("processForward")) {
      registry.add("forward/Tracks", "Number of forward track types", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    }
  }

  // This process functions fills the TracksWGTInBCs table.
  // It loops over all tracks. For the tracks with a 'good' timing it finds the associated BC.
  // If a track is ambiguous, then the associated BC is calculated with help of the trackTime.
  void processBarrel(BCs const& bcs, CCs const& collisions, TCs const& tracks, ATs const& ambTracks)
  {
    // run number
    if (bcs.size() <= 0) {
      return;
    }
    int rnum = bcs.iteratorAt(0).runNumber();

    // container to sort tracks with good timing according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> tracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all tracks and fill tracksInBCList
    LOGF(debug, "Number of barrel tracks: %d", tracks.size());
    for (auto const& track : tracks) {
      registry.get<TH1>(HIST("barrel/Tracks"))->Fill(0., 1.);

      // is this track a PV track?
      if (track.isPVContributor()) {
        registry.get<TH1>(HIST("barrel/Tracks"))->Fill(1., 1.);
      }

      // is this track an ambiguous track?
      auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
      if (ambTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("barrel/Tracks"))->Fill(3., 1.);
      } else {
        registry.get<TH1>(HIST("barrel/Tracks"))->Fill(2., 1.);
      }

      // only consider tracks with good timing
      if (std::abs(track.trackTimeRes()) <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("barrel/Tracks"))->Fill(4., 1.);

        // get first compatible BC
        if (ambTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("barrel/Tracks"))->Fill(5., 1.);

          // compute the BC closest in time
          auto firstCompatibleBC = ambTracksSlice.begin().bc().begin().globalBC();
          LOGF(debug, "First compatible BC %d Track time %f", firstCompatibleBC, track.trackTime());
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        } else {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          if (!track.has_collision()) {
            continue;
          }
          auto collision = track.collision_as<CCs>();
          if (!collision.has_foundBC()) {
            continue;
          }
          closestBC = collision.foundBC_as<BCs>().globalBC();
        }

        // update tracksInBCList
        LOGF(debug, "Updating tracksInBCList with %d", closestBC);
        tracksInBCList[closestBC].emplace_back((int32_t)track.globalIndex());
      }
    }

    // fill tracksWGTInBCs
    int indBCToStart = 0;
    int indBCToSave;
    for (auto const& tracksInBC : tracksInBCList) {
      LOGF(debug, "tracksInBC.first %d", tracksInBC.first);
      indBCToSave = -1;
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
      LOGF(debug, " BC %i/%u with %i tracks with good timing", indBCToSave, tracksInBC.first, tracksInBC.second.size());
      tracksWGTInBCs(indBCToSave, rnum, tracksInBC.first, tracksInBC.second);
    }
    LOGF(debug, "barrel done");
  }
  PROCESS_SWITCH(tracksWGTInBCs, processBarrel, "Process barrel tracks", false);

  // This process functions fills the FwdTracksWGTInBCs table.
  // It loops over all forward tracks. For the tracks with a 'good' timing it finds the associated BC.
  // If a track is ambiguous, then the associated BC is calculated with help of the trackTime.
  void processForward(BCs& bcs, CCs& collisions, aod::FwdTracks& fwdTracks, aod::AmbiguousFwdTracks& ambFwdTracks)
  {
    // run number
    if (bcs.size() <= 0) {
      return;
    }
    int rnum = bcs.iteratorAt(0).runNumber();

    // container to sort forward tracks according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> fwdTracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all forward tracks and fill fwdTracksInBCList
    LOGF(debug, "Number of forward tracks: %d", fwdTracks.size());
    for (auto const& fwdTrack : fwdTracks) {
      registry.get<TH1>(HIST("forward/Tracks"))->Fill(0., 1.);
      // is this track an ambiguous track?
      auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
      if (ambFwdTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("forward/Tracks"))->Fill(2., 1.);
      } else {
        registry.get<TH1>(HIST("forward/Tracks"))->Fill(1., 1.);
      }

      // only consider tracks with trackTimeRes < LHCBunchSpacingNS
      LOGF(debug, "Time resolution of fwdTrack %f", fwdTrack.trackTimeRes());
      if (std::abs(fwdTrack.trackTimeRes()) <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("forward/Tracks"))->Fill(3., 1.);

        // get first compatible BC
        // auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
        if (ambFwdTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("forward/Tracks"))->Fill(4., 1.);

          // compute the BC closest in time
          auto firstCompatibleBC = ambFwdTracksSlice.begin().bc().begin().globalBC();
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        } else {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          if (!fwdTrack.has_collision()) {
            continue;
          }
          auto collision = fwdTrack.collision_as<CCs>();
          if (!collision.has_foundBC()) {
            continue;
          }
          closestBC = collision.foundBC_as<BCs>().globalBC();
        }

        // update tracksInBCList
        LOGF(debug, "Updating fwdTracksInBCList with %d", closestBC);
        fwdTracksInBCList[closestBC].emplace_back((int32_t)fwdTrack.globalIndex());
      }
    }

    // fill fwdTracksWGTInBCs
    int indBCToStart = 0;
    int indBCToSave;
    for (auto const& fwdTracksInBC : fwdTracksInBCList) {
      indBCToSave = -1;
      // find corresponding BC
      for (auto ind = indBCToStart; ind < bcs.size(); ind++) {
        auto bc = bcs.rawIteratorAt(ind);
        if (bc.globalBC() == fwdTracksInBC.first) {
          indBCToSave = ind;
          indBCToStart = ind;
          break;
        }
        if (bc.globalBC() > fwdTracksInBC.first) {
          break;
        }
      }
      fwdTracksWGTInBCs(indBCToSave, rnum, fwdTracksInBC.first, fwdTracksInBC.second);
      LOGF(debug, " BC %i/%u with %i forward tracks with good timing", indBCToSave, fwdTracksInBC.first, fwdTracksInBC.second.size());
    }
    LOGF(debug, "forward done");
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
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDTracksLabels> outputTracksLabel;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  HistogramRegistry registry{
    "registry",
    {}};

  using TIBCs = aod::TracksWGTInBCs;
  using TIBC = TIBCs::iterator;
  using FTIBCs = aod::FwdTracksWGTInBCs;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, /*aod::TracksCov,*/ aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFbeta,
                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FTCs = aod::FwdTracks;

  SliceCache cache;
  Preslice<TCs> TCperCollision = aod::track::collisionId;
  Preslice<aod::FwdTracks> FWperCollision = aod::fwdtrack::collisionId;

  // update UDTables
  template <typename TTracks>
  void updateUDTables(bool onlyPV, int64_t colID, uint64_t bcnum, int rnum, float vx, float vy, float vz,
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
    outputCollsLabels(colID);

    // update DGTracks tables
    for (auto const& track : tracks) {
      if (track.isPVContributor() || !onlyPV) {
        updateUDTrackTables(outputCollisions.lastIndex(), track, bcnum);
      }
    }
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

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    if (context.mOptions.get<bool>("processTinBCs")) {
      registry.add("table/candCase", "#candCase", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    }

    if (context.mOptions.get<bool>("processFull")) {
      registry.add("data/candCase", "#candCase", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      registry.add("data/bcFlag", "#bcFlag", {HistType::kTH1F, {{64, -0.5, 63.5}}});
      registry.add("data/isDG1vsisDG2", "#isDG1vsisDG2", {HistType::kTH2F, {{13, -1.5, 11.5}, {13, -1.5, 11.5}}});
      registry.add("data/ntr1vsntr2All", "#ntr1vsntr2All", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("data/ntr1vsntr2Cand", "#ntr1vsntr2Cand", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
    }
  }

  // In this process function we run over all BCs in the TracksWGTInBCs table.
  // --> Only those BCs are considered!
  // If the BC has an associated collision, the information of this collision is extracted.
  // If the BC is not associated with a collision, then fill the UDtables with information availabl for the BC.
  void processTinBCs(TIBC const& tibc, BCs const& bcs, CCs const& collisions,
                     TCs const& tracks, aod::FwdTracks const& fwdtracks, FTIBCs const& ftibcs,
                     aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // fill FITInfo
    auto bcnum = tibc.bcnum();
    upchelpers::FITInfo fitInfo{};
    udhelpers::getFITinfo(fitInfo, bcnum, bcs, ft0s, fv0as, fdds);

    // check if DG event
    // distinguish between cases with and without associated BC
    // 1. candidate has associated BC and associated collision    -> vertex position: col.[posX(), posY(), posZ()]
    // 2. candidate has associated BC but no associated collision ->                  [-2., 2., -2.]
    // 3. candidate has no associated BC                          ->                  [-3., 3., -3.]
    int isDG = -1;
    float rtrwTOF = -1.;
    int8_t nCharge;
    if (tibc.has_bc()) {
      LOGF(debug, "[1.,2.] BC found");

      // get associated bc
      auto bc = tibc.bc_as<BCs>();

      // is there an associated collision?
      Partition<CCs> colSlize = aod::evsel::foundBCId == bc.globalIndex();
      colSlize.bindTable(collisions);

      if (colSlize.size() > 0) {
        LOGF(debug, "  1. BC has collision");
        colSlize.bindExternalIndices(&bcs);
        auto col = colSlize.begin();

        auto colTracks = tracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
        auto colFwdTracks = fwdtracks.sliceByCached(aod::fwdtrack::collisionId, col.globalIndex(), cache);
        auto bcRange = udhelpers::compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
        isDG = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);

        // update UDTables, case 1.
        if (isDG == 0) {
          registry.get<TH1>(HIST("table/candCase"))->Fill(1, 1.);
          rtrwTOF = udhelpers::rPVtrwTOF<true>(colTracks, col.numContrib());
          nCharge = udhelpers::netCharge<true>(colTracks);

          updateUDTables(false, col.globalIndex(), bc.globalBC(), bc.runNumber(), col.posX(), col.posY(), col.posZ(),
                         col.numContrib(), nCharge, rtrwTOF, colTracks, fitInfo);
        }
      } else {
        LOGF(debug, "  2. BC has NO collision");
        auto tracksArray = tibc.track_as<TCs>();
        auto bcRange = udhelpers::compatibleBCs(bc, bc.globalBC(), diffCuts.minNBCs(), bcs);

        // does BC have fwdTracks?
        if (ftibcs.size() > 0) {
          Partition<FTIBCs> ftibcSlice = aod::dgbcandidate::bcId == bc.globalIndex();
          ftibcSlice.bindTable(ftibcs);

          if (ftibcSlice.size() > 0) {
            ftibcSlice.bindExternalIndices(&fwdtracks);
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

        // update UDTables, case 2.
        if (isDG == 0) {
          registry.get<TH1>(HIST("table/candCase"))->Fill(2, 1.);
          rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
          nCharge = udhelpers::netCharge<false>(tracksArray);

          updateUDTables(false, -1, bc.globalBC(), bc.runNumber(), -2., 2., -2,
                         tracksArray.size(), nCharge, rtrwTOF, tracksArray, fitInfo);
        }
      }

      // fill UDZdcs
      if (isDG == 0 && bc.has_zdc()) {
        auto zdc = bc.zdc();
        std::vector<float> enes(zdc.energy()[0]);
        std::vector<uint8_t> chEs(zdc.channelE()[0]);
        std::vector<float> amps(zdc.amplitude()[0]);
        std::vector<float> times(zdc.time()[0]);
        std::vector<uint8_t> chTs(zdc.channelT()[0]);
        outputZdcs(outputCollisions.lastIndex(), enes, chEs, amps, times, chTs);
      }
    } else {
      LOGF(debug, "  3. BC NOT found");

      // the BC is not contained in the BCs table
      auto tracksArray = tibc.track_as<TCs>();
      auto bcRange = udhelpers::compatibleBCs(bcnum, diffCuts.minNBCs(), bcs);

      // does BC have fwdTracks?
      if (ftibcs.size() > 0) {
        Partition<FTIBCs> ftibcPart = aod::dgbcandidate::bcnum == bcnum;
        ftibcPart.bindTable(ftibcs);

        if (ftibcPart.size() > 0) {
          ftibcPart.bindExternalIndices(&fwdtracks);
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

      // update UDTables, case 3.
      if (isDG == 0) {
        registry.get<TH1>(HIST("table/candCase"))->Fill(3, 1.);
        rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
        nCharge = udhelpers::netCharge<false>(tracksArray);

        updateUDTables(false, -1, bcnum, tibc.runNumber(), -3., 3., -3,
                       tracksArray.size(), nCharge, rtrwTOF, tracksArray, fitInfo);
      }
    }
  }

  PROCESS_SWITCH(DGBCCandProducer, processTinBCs, "Produce UDTables", true);

  // In this process function runs over all BC mumbers
  // --> All possible BCs are considered!
  // 1. search for associated BCs and Collisions and extract the related information.
  // 2. search for entries in table TracksWGTInBCs and extract the related information.
  // ATTENTION: note that a BC can be included twice in the output tables. This can happen when
  //            the BC number is in both, the BCs and TracksWGTInBCs tables. The vertex position in
  //            UDCollisions is used to distinguish the cases.
  //              case:                                                           vertex position:
  //              BC in BCs table and has assoc. Collision:                       [posx(), posy(), posz()]
  //              BC in TracksWGTInBCs table:
  //                is also in BCs table and has assoc. Collision:                [-1., 1., -1]
  //                is also in BCs table and does not have assoc. Collision:      [-2., 2., -2.]
  //                is not in BCs table and hence does not have assoc. Collision: [-3., 3., -3.]
  void processFull(BCs const& bcs, CCs const& collisions,
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

    // run over all BC in bcs and tibcs
    int64_t lastCollision = 0;
    float vpos[3];
    upchelpers::FITInfo fitInfo{};
    auto col = collisions.iteratorAt(0);
    auto bc = bcs.iteratorAt(0);
    auto tibc = tibcs.iteratorAt(0);
    auto ftibc = ftibcs.iteratorAt(0);
    auto lastcol = collisions.iteratorAt(collisions.size() - 1);
    auto lastbc = bcs.iteratorAt(bcs.size() - 1);
    auto lasttibc = tibcs.iteratorAt(tibcs.size() - 1);
    auto lastftibc = ftibcs.iteratorAt(ftibcs.size() - 1);
    LOGF(debug, "collisions %d bcs %d tibcs %d ftibcs %d", collisions.size(), bcs.size(), tibcs.size(), ftibcs.size());

    // set first bcnum
    bool bc2go = bc != lastbc;
    bool tibc2go = tibcs.size() > 0 ? tibc != lasttibc : false;
    auto bcnum = bc.globalBC();
    if (tibc2go) {
      if (tibc.bcnum() < bcnum) {
        bcnum = tibc.bcnum();
      }
    }

    bool withCollision = false;
    while (bc2go || tibc2go) {
      LOGF(debug, "Testing bc %d/%d/%d", bcnum, bc.globalBC(), tibc.bcnum());
      // reset counters
      bcFlag = 1; // bit 0 is always set
      isDG1 = -1;
      isDG2 = -1;
      ntr1 = -1;
      ntr2 = -1;

      if (bc.globalBC() == bcnum) {
        SETBIT(bcFlag, 1);

        // find associated collision
        withCollision = false;
        auto goOn = col != lastcol;
        while (goOn) {
          if (col.has_foundBC()) {
            auto bc2u = col.foundBC_as<BCs>();
            if (bc2u.globalBC() >= bcnum) {
              goOn = false;
              if (bc2u.globalBC() == bcnum) {
                withCollision = true;
              }
            } else {
              col++;
            }
          } else {
            col++;
          }
          goOn &= col != lastcol;
        }
        LOGF(debug, "    withCollision %d", withCollision);

        if (withCollision) {
          // -> vertex position: col.[posX(), posY(), posZ()]
          SETBIT(bcFlag, 2);
          lastCollision = col.globalIndex();

          ntr1 = col.numContrib();
          auto bcRange = udhelpers::compatibleBCs(bc, bcnum, diffCuts.minNBCs(), bcs);
          auto colTracks = tracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
          auto colFwdTracks = fwdtracks.sliceByCached(aod::fwdtrack::collisionId, col.globalIndex(), cache);
          isDG1 = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);
          LOGF(debug, "  isDG1 %d with %d tracks", isDG1, ntr1);
          if (isDG1 == 0) {
            // this is a DG candidate with proper collision vertex
            SETBIT(bcFlag, 3);
            registry.get<TH1>(HIST("data/candCase"))->Fill(1, 1.);

            auto rtrwTOF = udhelpers::rPVtrwTOF<true>(colTracks, col.numContrib());
            auto nCharge = udhelpers::netCharge<true>(colTracks);
            udhelpers::getFITinfo(fitInfo, bcnum, bcs, ft0s, fv0as, fdds);
            updateUDTables(false, col.globalIndex(), bcnum, bc.runNumber(), col.posX(), col.posY(), col.posZ(),
                           col.numContrib(), nCharge, rtrwTOF, colTracks, fitInfo);
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
          }
        }
      }

      // find BC in TracksWGTInBCs table
      if (tibc2go) {
        if (tibc.bcnum() == bcnum) {
          SETBIT(bcFlag, 4);

          auto bcRange = udhelpers::compatibleBCs(bc, bcnum, diffCuts.minNBCs(), bcs);
          auto tracksArray = tibc.track_as<TCs>();
          ntr2 = tracksArray.size();

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

          LOGF(debug, "  isDG2 %d with %d tracks", isDG2, ntr2);
          if (isDG2 == 0) {
            // this is a DG candidate contained in TracksWGTInBCs
            SETBIT(bcFlag, 5);
            registry.get<TH1>(HIST("data/candCase"))->Fill(2, 1.);

            auto rtrwTOF = udhelpers::rPVtrwTOF<false>(tracksArray, tracksArray.size());
            auto nCharge = udhelpers::netCharge<false>(tracksArray);
            udhelpers::getFITinfo(fitInfo, bcnum, bcs, ft0s, fv0as, fdds);

            // distinguish different cases
            if (bc.globalBC() == bcnum) {
              if (withCollision) {
                vpos[0] = -1.;
                vpos[1] = 1.;
                vpos[2] = -1.;
              } else {
                vpos[0] = -2.;
                vpos[1] = 2.;
                vpos[2] = -2.;
              }
            } else {
              vpos[0] = -3.;
              vpos[1] = 3.;
              vpos[2] = -3.;
            }

            int64_t colID = withCollision ? col.globalIndex() : -1;
            updateUDTables(false, colID, bcnum, tibc.runNumber(), vpos[0], vpos[1], vpos[2],
                           tracksArray.size(), nCharge, rtrwTOF, tracksArray, fitInfo);
            // fill UDZdcs
            if (bc.globalBC() == bcnum) {
              if (bc.has_zdc()) {
                auto zdc = bc.zdc();
                auto enes = std::vector(zdc.energy().begin(), zdc.energy().end());
                auto chEs = std::vector(zdc.channelE().begin(), zdc.channelE().end());
                auto amps = std::vector(zdc.amplitude().begin(), zdc.amplitude().end());
                auto times = std::vector(zdc.time().begin(), zdc.time().end());
                auto chTs = std::vector(zdc.channelT().begin(), zdc.channelT().end());
                outputZdcs(outputCollisions.lastIndex(), enes, chEs, amps, times, chTs);
              }
            }
          }
        }
      }

      // update histograms
      registry.get<TH1>(HIST("data/bcFlag"))->Fill(bcFlag, 1.);
      registry.get<TH2>(HIST("data/isDG1vsisDG2"))->Fill(isDG1, isDG2);
      registry.get<TH2>(HIST("data/ntr1vsntr2All"))->Fill(ntr1, ntr2);
      if (isDG1 == 0 && isDG2 == 0) {
        registry.get<TH2>(HIST("data/ntr1vsntr2Cand"))->Fill(ntr1, ntr2);
      }

      // update bc and tibc
      if (bc2go) {
        if (bc.globalBC() == bcnum) {
          bc++;
        }
      }
      if (tibc2go) {
        if (tibc.bcnum() == bcnum) {
          tibc++;
        }
      }
      // determine next bcnum
      bc2go = bc != lastbc;
      tibc2go = tibcs.size() > 0 ? tibc != lasttibc : false;
      if (bc2go) {
        bcnum = bc.globalBC();
        if (tibc2go) {
          bcnum = std::min(bc.globalBC(), tibc.bcnum());
        }
      } else {
        if (tibc2go) {
          bcnum = tibc.bcnum();
        }
      }
    }
  }

  PROCESS_SWITCH(DGBCCandProducer, processFull, "Produce UDTables", true);

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
