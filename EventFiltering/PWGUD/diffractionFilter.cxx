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
// \brief A filter task for diffractive events
// \author P. Buehler, paul.buehler@oeaw.ac.at
// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/Core/DGCutparHolder.h"
#include "PWGUD/Core/DGSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Run 3
struct DGFilterRun3 {

  // Productions
  Produces<aod::DiffractionFilters> filterTable;

  // DGCutparHolders
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> diffCutsHolder{"DiffCuts", {}, "Diffractive events cuts"};

  // DG selector
  DGSelector dgSelector;

  // histograms with cut statistics
  // bin:
  //   1: All collisions
  //   2: DG candidate
  //   3: not clean FIT
  //   4: number of FwdTracks > 0
  //   5: not all global tracks are vtx tracks
  //   6: not all vtx tracks are global tracks
  //   7: fraction of tracks with TOF hit too low
  //   8: number of vtx tracks out of range
  //   9: has not good PID information
  //  10: track pt out of range
  //  11: track eta out of range
  //  12: net charge out of range
  //  13: IVM out of range
  HistogramRegistry registry{
    "registry",
    {
      {"aftercuts", "#aftercuts", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
    }};

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)diffCutsHolder;
  }

  // some general Collisions and Tracks filter
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  // using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
  //                       aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
  //                       aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  // using MFs = aod::MFTTracks;
  using FWs = aod::FwdTracks;

  void process(CC const& collision,
               BCs const& bcs,
               TCs& tracks,
               // MFs& mfttracks,
               FWs& fwdtracks,
               aod::Zdcs& zdcs,
               aod::FT0s& ft0s,
               aod::FV0As& fv0as,
               aod::FDDs& fdds)
  {

    // initialize
    LOGF(debug, "<DGFilterRun3. Collision %d", collision.globalIndex());
    bool ccs{false};
    registry.fill(HIST("aftercuts"), 0.);

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
    LOGF(debug, "  Number of compatible BCs in +- %i / %i dtcoll: %i", diffCuts.NDtcoll(), diffCuts.minNBCs(), bcRange.size());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, tracks, fwdtracks);

    // update after cut histogram
    registry.fill(HIST("aftercuts"), isDGEvent + 1);

    // update filterTable
    ccs = (isDGEvent == 0);
    filterTable(ccs);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGFilterRun3>(cfgc, TaskName{"DGfilterRun3"}),
  };
}
