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
// O2 includes
//
// \brief A filter task for diffractive BCs
// \author P. Buehler, paul.buehler@oeaw.ac.at
// \since December, 2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/TableProducer/DGBCCandProducer.h"
#include "PWGUD/Core/DGCutparHolder.h"
#include "PWGUD/Core/DGSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "../filterTables.h"

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
    uint64_t closestBC;

    // loop over all tracks and fill tracksInBCList
    LOGF(debug, "Number of barrel tracks: %d", tracks.size());
    for (auto const& track : tracks) {
      registry.get<TH1>(HIST("barrelTracks"))->Fill(0., 1.);
      auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
      if (ambTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(2., 1.);
      } else {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(1., 1.);
      }

      // only consider tracks with good timing
      LOGF(debug, "Track time %f resolution %f", track.trackTime(), track.trackTimeRes());
      if (track.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("barrelTracks"))->Fill(3., 1.);

        // get first compatible BC
        // auto ambTracksSlice = ambTracks.sliceBy(perTrack, track.globalIndex());
        LOGF(debug, "Size of ambTracksSlice %d", ambTracksSlice.size());
        if (ambTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("barrelTracks"))->Fill(4., 1.);

          // compute the BC closest in time
          auto atr = ambTracksSlice.begin();
          LOGF(debug, "  Number of BCs %d", atr.bc().size());
          if (atr.bc().size() > 0) {
            auto firstCompatibleBC = atr.bc().begin().globalBC();
            closestBC = (uint64_t)(firstCompatibleBC +
                                   (track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));
            // update tracksInBCList
            LOGF(debug, "closestBC %d", closestBC);
            tracksInBCList[closestBC].emplace_back((int32_t)track.globalIndex());
          }
        } else {
          // registry.get<TH1>(HIST("barrelTracks"))->Fill(2., 1.);

          // this track is not ambiguous, has hence a unique association to a collision/BC
          auto col = track.collision_as<CCs>();
          LOGF(debug, "  has BC %d", col.has_foundBC());
          closestBC = track.collision_as<CCs>().foundBC_as<BCs>().globalBC();

          // update tracksInBCList
          LOGF(debug, "closestBC %d", closestBC);
          tracksInBCList[closestBC].emplace_back((int32_t)track.globalIndex());
        }
      }
      LOGF(debug, "track finished.\n");
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
  PROCESS_SWITCH(tracksWGTInBCs, processBarrel, "Process barrel tracks", true);

  void processForward(BCs& bcs, CCs& collisions, aod::FwdTracks& fwdTracks, aod::AmbiguousFwdTracks& ambFwdTracks)
  {
    // run number
    int rnum = bcs.iteratorAt(0).runNumber();

    // container to sort forward tracks according to their matching/closest BC
    std::map<uint64_t, std::vector<int32_t>> fwdTracksWGTInBCList{};
    uint64_t closestBC = 0;

    // loop over all forward tracks and fill fwdTracksWGTInBCList
    LOGF(debug, "Number of forward tracks: %d", fwdTracks.size());
    for (auto const& fwdTrack : fwdTracks) {
      registry.get<TH1>(HIST("forwardTracks"))->Fill(0., 1.);
      auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
      if (ambFwdTracksSlice.size() > 0) {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(2., 1.);
      } else {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(1., 1.);
      }

      // only consider tracks with trackTimeRes < LHCBunchSpacingNS
      LOGF(debug, "FwdTrack time %f resolution %f", fwdTrack.trackTime(), fwdTrack.trackTimeRes());
      if (fwdTrack.trackTimeRes() <= o2::constants::lhc::LHCBunchSpacingNS) {
        registry.get<TH1>(HIST("forwardTracks"))->Fill(3., 1.);

        // get first compatible BC
        // auto ambFwdTracksSlice = ambFwdTracks.sliceBy(perFwdTrack, fwdTrack.globalIndex());
        LOGF(debug, "Size of ambFwdTracksSlice %d", ambFwdTracksSlice.size());
        if (ambFwdTracksSlice.size() > 0) {
          registry.get<TH1>(HIST("forwardTracks"))->Fill(4., 1.);

          // compute the BC closest in time
          auto aftr = ambFwdTracksSlice.begin();
          LOGF(debug, "  Number of BCs %d", aftr.bc().size());
          if (aftr.bc().size() > 0) {
            auto firstCompatibleBC = aftr.bc().begin().globalBC();
            closestBC = (uint64_t)(firstCompatibleBC +
                                   (fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS));

            // update fwdTracksWGTInBCList
            LOGF(debug, "closestBC %d", closestBC);
            fwdTracksWGTInBCList[closestBC].emplace_back((int32_t)fwdTrack.globalIndex());
          }
        } else {
          // registry.get<TH1>(HIST("forwardTracks"))->Fill(2., 1.);

          // this track is not ambiguous, has hence a unique association to a collision/BC
          auto col = fwdTrack.collision_as<CCs>();
          LOGF(debug, "  has BC %d", col.has_foundBC());
          closestBC = col.foundBC_as<BCs>().globalBC();

          // update fwdTracksWGTInBCList
          LOGF(debug, "closestBC %d", closestBC);
          fwdTracksWGTInBCList[closestBC].emplace_back((int32_t)fwdTrack.globalIndex());
        }
      }
      LOGF(debug, "FwdTrack finished.\n");
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
  PROCESS_SWITCH(tracksWGTInBCs, processForward, "Process forward tracks", true);
};

// -----------------------------------------------------------------------------
// Run 3
struct DGBCFilterRun3 {

  // Productions
  Produces<aod::DiffractionBCFilters> bcfilterTable;

  // DGCutparHolders
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  // histograms with cut statistics
  // bin:
  //   1: All collisions
  //   2: DG candidate
  //   3: FIT veto
  //   4: number of FwdTracks > 0
  //   5: not all global tracks are PV tracks
  //   6: not all PV tracks are global tracks
  //   7: ITS only PV tracks
  //   8: fraction of tracks with TOF hit too low
  //   9: number of PV tracks out of range
  //  10: PV tracks without good PID information
  //  11: PV track pt out of range
  //  13: PV track eta out of range
  //  14: net charge out of range
  //  15: IVM out of range
  HistogramRegistry registry{
    "registry",
    {
      {"aftercut", "#aftercut", {HistType::kTH1F, {{14, -0.5, 13.5}}}},
    }};

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
  }

  // some general Collisions and Tracks filter
  using TIBCs = aod::TracksWGTInBCs;
  using TIBC = TIBCs::iterator;
  using FTIBCs = aod::FwdTracksWGTInBCs;
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FTCs = aod::FwdTracks;

  // using MFs = aod::MFTTracks;
  using FWs = aod::FwdTracks;

  void process(BCs const& bcs,
               CCs const& collisions,
               TCs& tracks,
               // MFs& mfttracks,
               FWs& fwdtracks,
               TIBCs const& tibcs, FTIBCs const& ftibcs,
               aod::Zdcs& zdcs,
               aod::FT0s& ft0s,
               aod::FV0As& fv0as,
               aod::FDDs& fdds)
  {

    // Advance these pointers step-by-step
    auto tibc = tibcs.iteratorAt(0);
    auto ftibc = ftibcs.iteratorAt(0);
    auto lasttibc = tibcs.iteratorAt(tibcs.size() - 1);
    auto lastftibc = ftibcs.iteratorAt(ftibcs.size() - 1);

    // loop over bcs
    int isDGBC;
    for (auto bc : bcs) {
      auto bcnum = bc.globalBC();
      auto ccs = false;

      // find BC in TIBCs table
      while (tibc.bcnum() < bcnum && tibc != lasttibc) {
        ++tibc;
      }
      if (tibc.bcnum() == bcnum) {
        registry.fill(HIST("aftercut"), 0.);

        // get tracks associated with BC
        auto tracksArray = tibc.track_as<TCs>();

        // obtain slice of compatible BCs
        auto bcRange = udhelpers::compatibleBCs(bc, bcnum, diffCuts.minNBCs(), bcs);

        // find BC in FTIBCs table
        while (ftibc.bcnum() < bcnum && ftibc != lastftibc) {
          ++ftibc;
        }

        // apply DG selection
        if (ftibc.bcnum() == bcnum) {
          auto fwdTracksArray = ftibc.fwdtrack_as<FTCs>();
          isDGBC = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        } else {
          auto fwdTracksArray = FTCs{{fwdtracks.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
          isDGBC = dgSelector.IsSelected(diffCuts, bcRange, tracksArray, fwdTracksArray);
        }

        // save decision
        ccs = (isDGBC == 0);
        if (ccs) {
          LOGF(debug, "BC %d is a DG candidate!", bcnum);
        }

        // update after cut histogram
        registry.fill(HIST("aftercut"), isDGBC + 1);
      }

      // update filterTable
      bcfilterTable(ccs);
    }
  }
};

// -----------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tracksWGTInBCs>(cfgc, TaskName{"trackswgtinbcs"}),
    adaptAnalysisTask<DGBCFilterRun3>(cfgc, TaskName{"DGBCfilterRun3"}),
  };
}

// -----------------------------------------------------------------------------
