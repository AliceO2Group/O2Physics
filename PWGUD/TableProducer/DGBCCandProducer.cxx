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
#include "PWGUD/Core/UDHelperFunctions.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"
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
  }

  void processBarrel(BCs const& bcs, CCs const& collisions, TCs const& tracks, ATs const& ambTracks)
  {
    // container to sort tracks with good timing according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> tracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all tracks and fill tracksInBCList
    for (auto const& track : tracks) {
      // only consider tracks with good timing
      if (track.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {

        // get first compatible BC
        auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
        if (ambTracksSlice.size() > 0) {
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
        tracksWGTInBCs(indBCToSave, tracksInBC.first, tracksInBC.second);
        LOGF(debug, " BC %i/%u with %i tracks with good timing", indBCToSave, tracksInBC.first, tracksInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWGTInBCs, processBarrel, "Process barrel tracks", true);

  void processForward(BCs& bcs, CCs& collisions, aod::FwdTracks& fwdTracks, aod::AmbiguousFwdTracks& ambFwdTracks)
  {
    // container to sort forward tracks according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> fwdTracksWGTInBCList{};
    uint64_t closestBC = 0;

    // loop over all forward tracks and fill fwdTracksWGTInBCList
    for (auto const& fwdTrack : fwdTracks) {
      // only consider tracks with trackTimeRes < LHCBunchSpacingNS
      LOGF(info, "Time resolution of fwdTrack %f", fwdTrack.trackTimeRes());
      if (fwdTrack.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {

        // get first compatible BC
        auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
        if (ambFwdTracksSlice.size() > 0) {
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
        fwdTracksWGTInBCs(indBCToSave, fwdTracksWGTInBC.first, fwdTracksWGTInBC.second);
        LOGF(debug, " BC %i/%u with %i forward tracks with good timing", indBCToSave, fwdTracksWGTInBC.first, fwdTracksWGTInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWGTInBCs, processForward, "Process forward tracks", true);
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
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FTCs = aod::FwdTracks;

  Preslice<CCs> CCperBC = aod::evsel::foundBCId;
  Preslice<TCs> TCperCollision = aod::track::collisionId;
  Preslice<aod::FwdTracks> FWperCollision = aod::fwdtrack::collisionId;
  Preslice<TIBCs> TIBCperBC = aod::dgbcandidate::bcId;
  Preslice<FTIBCs> FTIBCperBC = aod::dgbcandidate::bcId;

  // extract FIT information
  upchelpers::FITInfo getFITinfo(BC const& bc, BCs const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // FITinfo
    upchelpers::FITInfo info{};

    // FT0
    if (bc.has_foundFT0()) {
      auto ft0 = bc.foundFT0();
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
      auto fv0a = bc.foundFV0();
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
      auto fdd = bc.foundFDD();
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

    // fill BG and BB flags in adjacent BCs [-15, 15]
    auto minbcInd = bc.globalIndex() > 15 ? bc.globalIndex() - 15 : 0;
    auto maxbcInd = bc.globalIndex() < bcs.size() - 1 - 15 ? bc.globalBC() + 15 : bcs.size() - 1;
    auto minbc = bc.globalBC() - 15;
    auto maxbc = bc.globalBC() + 15;
    for (uint ind = minbcInd; ind <= maxbcInd; ind++) {

      // is this a relevant BC
      auto bc2u = bcs.iteratorAt(ind);
      if (bc2u.globalBC() > maxbc) {
        break;
      }
      if (bc2u.globalBC() < minbc) {
        continue;
      }

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

    return info;
  }

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack, typename TBC>
  void updateUDTrackTables(TTrack const& track, TBC const& bc)
  {
    outputTracks(outputCollisions.lastIndex(), track.px(), track.py(), track.pz(), track.sign(),
                 bc.globalBC(), track.trackTime(), track.trackTimeRes());
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
    outputTracksFlag(track.has_collision() ? track.collisionId() : -1,
                     track.isPVContributor());
  }

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    if (context.mOptions.get<bool>("processQA")) {
      registry.add("isDG1vsisDG2", "#isDG1vsisDG2", {HistType::kTH2F, {{13, -1.5, 11.5}, {13, -1.5, 11.5}}});
      registry.add("ntr1vsntr2All", "#ntr1vsntr2All", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("ntr1vsntr2Cand", "#ntr1vsntr2Cand", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("ptvsdcaxy", "#ptvsdcaxy", {HistType::kTH2F, {{50, 0., 5.}, {50, -5., 5.}}});
    }
  }

  void process(TIBC const& tibc, BCs const& bcs, CCs const& collisions,
               TCs const& tracks, aod::FwdTracks const& fwdtracks, FTIBCs const& ftibcs,
               aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {

    // leave if tibc has no associated BC
    if (!tibc.has_bc()) {
      return;
    }

    // get bc
    auto bc = tibc.bc_as<BCs>();

    // fill FITInfo
    upchelpers::FITInfo fitInfo = getFITinfo(bc, bcs, ft0s, fv0as, fdds);

    // check if DG event
    float rtrwTOF = -1.;
    int8_t nCharge;

    // is there an associated collision?
    int isDG = -1;
    auto colSlize = collisions.sliceBy(CCperBC, bc.globalIndex());
    if (colSlize.size() > 0) {
      auto col = colSlize.rawIteratorAt(0);
      auto colTracks = tracks.sliceBy(TCperCollision, col.globalIndex());
      auto colFwdTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
      auto bcRange = compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
      isDG = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);

      // update UDTables
      if (isDG == 0) {
        rtrwTOF = rPVtrwTOF<true>(colTracks, col.numContrib());
        nCharge = netCharge(colTracks);

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
          if (track.isPVContributor()) {
            updateUDTrackTables(track, bc);
          }
        }
      }
    } else {
      auto tracksArray = tibc.track_as<TCs>();
      auto bcRange = compatibleBCs(bc, bc.globalBC(), diffCuts.minNBCs(), bcs);

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
        rtrwTOF = rPVtrwTOF<false>(tracksArray, tracksArray.size());
        nCharge = netCharge(tracksArray);

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
          updateUDTrackTables(track, bc);
        }
      }
    }
  }

  void processQA(BCs const& bcs, CCs const& collisions,
                 TCs const& tracks, FTCs const& fwdtracks, TIBCs const& tibcs, FTIBCs const& ftibcs,
                 aod::Zdcs const& zdcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {

    // loop over BCs
    int isDG1, isDG2;
    int ntr1, ntr2;
    for (auto bc : bcs) {
      // reset counters
      isDG1 = -1;
      isDG2 = -1;
      ntr1 = -1;
      ntr2 = -1;

      // check for associated collision
      auto colSlize = collisions.sliceBy(CCperBC, bc.globalIndex());
      if (colSlize.size() > 0) {
        auto col = colSlize.begin();
        auto colTracks = tracks.sliceBy(TCperCollision, col.globalIndex());
        auto colFwdTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
        auto bcRange = compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
        isDG1 = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFwdTracks);
        ntr1 = col.numContrib();
      }

      // find corresponding entry in tibcs
      auto tibcSlice = tibcs.sliceBy(TIBCperBC, bc.globalIndex());
      if (tibcSlice.size() > 0) {
        auto tibc = tibcSlice.begin();

        // check collision to be DGCandidate -> isDG2
        auto bcRange = compatibleBCs(bc, bc.globalBC(), diffCuts.minNBCs(), bcs);
        auto tracksArray = tibc.track_as<TCs>();
        auto ftibcSlice = ftibcs.sliceBy(FTIBCperBC, bc.globalIndex());
        if (ftibcSlice.size() > 0) {
          auto fwdTracksArray = ftibcSlice.begin().fwdtrack_as<FTCs>();
          isDG2 = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        } else {
          auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
          isDG2 = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        }
        ntr2 = tracksArray.size();

        // update ptvsdcaxy
        for (auto const& track : tracksArray) {
          registry.get<TH2>(HIST("ptvsdcaxy"))->Fill(track.pt(), track.dcaXY());
        }
      }

      // update histogram isDG1vsisDG2 and ntr1vsntr2xxx
      registry.get<TH2>(HIST("isDG1vsisDG2"))->Fill(isDG1, isDG2);
      registry.get<TH2>(HIST("ntr1vsntr2All"))->Fill(ntr1, ntr2);
      if (isDG1 == 0 && isDG2 == 0) {
        registry.get<TH2>(HIST("ntr1vsntr2Cand"))->Fill(ntr1, ntr2);
      }
    }
  }

  PROCESS_SWITCH(DGBCCandProducer, processQA, "Produce QA histograms for DGBCCandProducer", false);
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
