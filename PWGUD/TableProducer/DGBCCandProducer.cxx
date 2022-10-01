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

  void init(InitContext& context)
  {
  }

  // ---------------------------------------------------------------------------

  Preslice<aod::AmbiguousTracks> perTrack = aod::ambiguous::trackId;

  void process(TCs& tracks, BCs& bcs, CCs& collisions, ATs& ambTracks)
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
                                 (track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + 0.5));
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
};

// -----------------------------------------------------------------------------
struct DGBCCandProducer {
  // data tables
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTrackCollisionIDs> outputTracksCollisionsId;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext& context)
  {
    diffCuts = (DGCutparHolder)DGCuts;

    if (context.mOptions.get<bool>("processQA")) {
      registry.add("isDG1vsisDG2", "#isDG1vsisDG2", {HistType::kTH2F, {{13, -1.5, 11.5}, {13, -1.5, 11.5}}});
      registry.add("ntr1vsntr2", "#ntr1vsntr2", {HistType::kTH2F, {{52, -1.5, 50.5}, {52, -1.5, 50.5}}});
    }
  }

  using TWTs = aod::TrackswTOFInBCs;
  using TWT = TWTs::iterator;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  Preslice<CCs> CCperBC = aod::evsel::foundBCId;
  Preslice<TCs> TCperCollision = aod::track::collisionId;
  Preslice<FWs> FWperCollision = aod::fwdtrack::collisionId;
  Preslice<TWTs> TWTperBC = aod::dgbcandidate::bcId;

  // function to update UDTracks, UDTracksPID, UDTracksExtra, and UDTrackCollisionIDs
  template <typename TTrack, typename TBC>
  void updateUDTrackTables(TTrack const& track, TBC const& bc)
  {
    outputTracks(track.px(), track.py(), track.pz(), track.sign(),
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
    outputTracksCollisionsId(outputCollisions.lastIndex());
  }

  void process(TWT& twt, BCs& bcs, CCs& collisions, TCs& tracks, FWs& fwdtracks,
               aod::Zdcs& zdcs, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {

    // leave if twt has no associated BC
    if (!twt.has_bc()) {
      return;
    }

    // get bc and related tracks
    auto bc = twt.bc_as<BCs>();

    // check if DG event
    float rtrwTOF = -1.;
    int8_t nCharge;

    // is there an associated collision?
    auto colSlize = collisions.sliceBy(CCperBC, bc.globalIndex());
    if (colSlize.size() > 0) {
      auto col = colSlize.rawIteratorAt(0);
      auto colTracks = tracks.sliceBy(TCperCollision, col.globalIndex());
      auto colFWDTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
      auto bcRange = compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
      if (dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFWDTracks) == 0) {
        rtrwTOF = rPVtrwTOF(colTracks, col.numContrib());
        nCharge = netCharge(colTracks);

        outputCollisions(bc.globalBC(), bc.runNumber(),
                         col.posX(), col.posY(), col.posZ(),
                         col.numContrib(), nCharge,
                         rtrwTOF,
                         0., 0., 0., 0., 0);

        // update DGTracks tables
        for (auto& track : colTracks) {
          if (track.isPVContributor()) {
            updateUDTrackTables(track, bc);
          }
        }
      }
    } else {
      auto tracksArray = twt.track_as<TCs>();
      if (dgSelector.IsSelected(diffCuts, bc, tracksArray) == 0) {
        rtrwTOF = -1.;
        nCharge = netCharge(tracksArray);

        outputCollisions(bc.globalBC(), bc.runNumber(),
                         0., 0., 0.,
                         tracksArray.size(), nCharge,
                         rtrwTOF,
                         0., 0., 0., 0., 0);

        // update DGTracks tables
        for (auto& track : tracksArray) {
          updateUDTrackTables(track, bc);
        }
      }
    }
  }

  void processQA(TWTs& twts, BCs& bcs, CCs& collisions, TCs& tracks, FWs& fwdtracks,
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
        auto colFWDTracks = fwdtracks.sliceBy(FWperCollision, col.globalIndex());
        auto bcRange = compatibleBCs(col, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
        isDG1 = dgSelector.IsSelected(diffCuts, col, bcRange, colTracks, colFWDTracks);
        ntr1 = col.numContrib();
      }

      // find corresponding entry in twts
      auto twtSlice = twts.sliceBy(TWTperBC, bc.globalIndex());
      if (twtSlice.size() > 0) {
        auto twt = twtSlice.begin();

        // check collision to be DGCandidate -> isDG2
        auto tracksArray = twt.track_as<TCs>();
        isDG2 = dgSelector.IsSelected(diffCuts, bc, tracksArray);
        ntr2 = tracksArray.size();
      }

      // update histogram isDG1vsisDG2 and ntr1vsntr2
      registry.get<TH2>(HIST("isDG1vsisDG2"))->Fill(isDG1, isDG2);
      if (isDG1 == 0 && isDG2 == 0) {
        registry.get<TH2>(HIST("ntr1vsntr2"))->Fill(ntr1, ntr2);
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
