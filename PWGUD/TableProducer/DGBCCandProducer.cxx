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
#include "PWGUD/DataModel/UDTables.h"
#include "DGBCCandProducer.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// -----------------------------------------------------------------------------
struct tracksWTOFInBCs {
  Produces<aod::TrackswTOFInBCs> tracksWTOFInBCs;
  Produces<aod::FwdTracksInBCs> fwdTracksInBCs;

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

  void processBarrel(TCs& tracks, BCs& bcs, CCs& collisions, ATs& ambTracks)
  {
    // container to sort tracks with TOF hit according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> tracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all tracks and fill tracksInBCList
    for (auto& track : tracks) {
      // only consider tracks with TOF hit
      if (track.hasTOF()) {

        // get first compatible BC
        auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
        if (ambTracksSlice.size() == 0) {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          closestBC = track.collision_as<CCs>().bc_as<BCs>().globalBC();
        } else {
          // compute the BC closest in time
          auto firstCompatibleBC = ambTracksSlice.begin().bc().begin().globalBC();
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        }

        // update tracksInBCList
        tracksInBCList[closestBC].emplace_back((int32_t)track.globalIndex());
      }
    }

    // fill tracksWTOFInBCs
    int indBCToSave;
    for (auto tracksInBC : tracksInBCList) {
      indBCToSave = -1;
      if (tracksInBC.second.size() > 0) {
        // find corresponding BC
        for (auto ind = 0; ind < bcs.size(); ind++) {
          auto bc = bcs.rawIteratorAt(ind);
          if (bc.globalBC() == tracksInBC.first) {
            indBCToSave = ind;
            break;
          }
          if (bc.globalBC() > tracksInBC.first) {
            break;
          }
        }
        tracksWTOFInBCs(indBCToSave, tracksInBC.first, tracksInBC.second);
        LOGF(debug, " BC %i/%u with %i tracks with TOF", indBCToSave, tracksInBC.first, tracksInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWTOFInBCs, processBarrel, "Process barrel tracks", true);

  void processForward(aod::FwdTracks& fwdTracks, BCs& bcs, CCs& collisions, aod::AmbiguousFwdTracks& ambFwdTracks)
  {
    // container to sort forward tracks according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> fwdTracksInBCList{};
    uint64_t closestBC = 0;

    // loop over all forward tracks and fill fwdTracksInBCList
    for (auto& fwdTrack : fwdTracks) {
      // only consider tracks with trackTimeRes < LHCBunchSpacingNS
      LOGF(debug, "fwdTrack.trackTimeRes %f", fwdTrack.trackTimeRes());
      if (fwdTrack.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {

        // get first compatible BC
        auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
        if (ambFwdTracksSlice.size() == 0) {
          // this track is not ambiguous, has hence a unique association to a collision/BC
          closestBC = fwdTrack.collision_as<CCs>().bc_as<BCs>().globalBC();
        } else {
          // compute the BC closest in time
          auto firstCompatibleBC = ambFwdTracksSlice.begin().bc().begin().globalBC();
          closestBC = (uint64_t)(firstCompatibleBC +
                                 (fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
        }

        // update tracksInBCList
        fwdTracksInBCList[closestBC].emplace_back((int32_t)fwdTrack.globalIndex());
      }
    }

    // fill fwdTracksInBCs
    int indBCToSave;
    for (auto fwdTracksInBC : fwdTracksInBCList) {
      indBCToSave = -1;
      if (fwdTracksInBC.second.size() > 0) {
        // find corresponding BC
        for (auto ind = 0; ind < bcs.size(); ind++) {
          auto bc = bcs.rawIteratorAt(ind);
          if (bc.globalBC() == fwdTracksInBC.first) {
            indBCToSave = ind;
            break;
          }
          if (bc.globalBC() > fwdTracksInBC.first) {
            break;
          }
        }
        fwdTracksInBCs(indBCToSave, fwdTracksInBC.first, fwdTracksInBC.second);
        LOGF(debug, " BC %i/%u with %i tracks with TOF", indBCToSave, fwdTracksInBC.first, fwdTracksInBC.second.size());
      }
    }
  }
  PROCESS_SWITCH(tracksWTOFInBCs, processForward, "Process forward tracks", true);
};

// -----------------------------------------------------------------------------
struct DGBCCandProducer {
  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  HistogramRegistry registry{
    "registry",
    {}};

  using TIBCs = aod::TrackswTOFInBCs;
  using TIBC = TIBCs::iterator;
  using FTIBCs = aod::FwdTracksInBCs;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FTCs = aod::FwdTracks;

  Preslice<CCs> CCperBC = aod::evsel::foundBCId;
  Preslice<TCs> TCperCollision = aod::track::collisionId;
  Preslice<FTCs> FWperCollision = aod::fwdtrack::collisionId;
  Preslice<TIBCs> TIBCperBC = aod::dgbcandidate::bcId;
  Preslice<FTIBCs> FTIBCperBC = aod::dgbcandidate::bcId;

  // function to update UDTracks, UDTracksPID, UDTracksExtra, and UDTrackCollisionIDs
  template <typename TTrack, typename TBC>
  void updateUDTrackTables(TTrack const& track, TBC const& bc)
  {
    outputTracks(outputCollisions.lastIndex(), track.px(), track.py(), track.pz(), track.sign(),
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
  }

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    if (context.mOptions.get<bool>("processQA")) {
      registry.add("isDG1vsisDG2", "#isDG1vsisDG2", {HistType::kTH2F, {{13, -1.5, 11.5}, {13, -1.5, 11.5}}});
      registry.add("ntr1vsntr2All", "#ntr1vsntr2All", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
      registry.add("ntr1vsntr2Cand", "#ntr1vsntr2Cand", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
    }
  }

  void process(TIBC& tibc, FTIBCs& ftibcs, BCs& bcs, CCs& collisions, TCs& tracks, FTCs& fwdtracks,
               aod::Zdcs& zdcs, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {

    // leave if tibc has no associated BC
    if (!tibc.has_bc()) {
      return;
    }

    // get bc
    auto bc = tibc.bc_as<BCs>();

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
        rtrwTOF = rPVtrwTOF(colTracks, col.numContrib());
        nCharge = netCharge(colTracks);

        outputCollisions(bc.globalBC(), bc.runNumber(),
                         col.posX(), col.posY(), col.posZ(),
                         col.numContrib(), nCharge,
                         rtrwTOF);

        // update DGTracks tables
        for (auto& track : colTracks) {
          if (track.isPVContributor()) {
            updateUDTrackTables(track, bc);
          }
        }
      }
    } else {
      auto tracksArray = tibc.track_as<TCs>();

      // does BC have fwdTracks?
      auto ftibcSlice = ftibcs.sliceBy(FTIBCperBC, bc.globalIndex());
      if (ftibcSlice.size() > 0) {
        auto fwdTracksArray = ftibcSlice.begin().fwdtrack_as<FTCs>();
        isDG = dgSelector.IsSelected(diffCuts, bc, tracksArray, fwdTracksArray);
      } else {
        auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
        isDG = dgSelector.IsSelected(diffCuts, bc, tracksArray, fwdTracksArray);
      }

      // update UDTables
      if (isDG == 0) {
        rtrwTOF = -1.;
        nCharge = netCharge(tracksArray);

        outputCollisions(bc.globalBC(), bc.runNumber(),
                         0., 0., 0.,
                         tracksArray.size(), nCharge,
                         rtrwTOF);

        // update DGTracks tables
        for (auto& track : tracksArray) {
          updateUDTrackTables(track, bc);
        }
      }
    }
  }

  void processQA(TIBCs& tibcs, FTIBCs& ftibcs, BCs& bcs, CCs& collisions, TCs& tracks, FTCs& fwdtracks,
                 aod::Zdcs& zdcs, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
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
        auto tracksArray = tibc.track_as<TCs>();
        auto ftibcSlice = ftibcs.sliceBy(FTIBCperBC, bc.globalIndex());
        if (ftibcSlice.size() > 0) {
          auto fwdTracksArray = ftibcSlice.begin().fwdtrack_as<FTCs>();
          isDG2 = dgSelector.IsSelected(diffCuts, bc, tracksArray, fwdTracksArray);
        } else {
          auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
          isDG2 = dgSelector.IsSelected(diffCuts, bc, tracksArray, fwdTracksArray);
        }
        ntr2 = tracksArray.size();
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
    adaptAnalysisTask<tracksWTOFInBCs>(cfgc, TaskName{"trackswtofinbcs"}),
    adaptAnalysisTask<DGBCCandProducer>(cfgc, TaskName{"dgbccandproducer"}),
  };
}

// -----------------------------------------------------------------------------
