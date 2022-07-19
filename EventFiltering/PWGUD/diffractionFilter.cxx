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
//
//        options:  X in [2pi, 4pi, 2K, 4K]
//
//               DiffCutsX.mNDtcoll(4)
//               DiffCutsX.mMinNTracks(0)
//               DiffCutsX.mMaxNTracks(10000)
//               DiffCutsX.mMinNetCharge(0)
//               DiffCutsX.mMaxNetCharge(0)
//               DiffCutsX.mPidHypo(211)
//               DiffCutsX.mMinPosz(-1000.)
//               DiffCutsX.mMaxPosz(1000.)
//               DiffCutsX.mMinPt(0.)
//               DiffCutsX.mMaxPt(1000.)
//               DiffCutsX.mMinEta(-1.)
//               DiffCutsX.mMaxEta(1.)
//               DiffCutsX.mMinIVM(0.)
//               DiffCutsX.mMaxIVM(1000.)
//               DiffCutsX.mMaxnSigmaTPC(1000.)
//               DiffCutsX.mMaxnSigmaTOF(1000.)
//               DiffCutsX.mFITAmpLimits({0., 0., 0., 0., 0.})
//
//        usage: copts="--configuration json://DiffFilterConfig.json -b"
//               kopts="--aod-writer-keep dangling --aod-writer-resfile DiffSelection"
//
//               o2-analysis-timestamp $copts |
//               o2-analysis-track-propagation $copts |
//               o2-analysis-event-selection $copts |
//               o2-analysis-multiplicity-table $copts |
//               o2-analysis-trackextension $copts |
//               o2-analysis-trackselection $copts |
//               o2-analysis-pid-tpc-full $copts |
//               o2-analysis-pid-tof-base $copts |
//               o2-analysis-pid-tof-full $copts |
//               o2-analysis-diffraction-filter $copts $kopts > diffractionFilter.log

// \author P. Buehler, paul.buehler@oeaw.ac.at
// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "DGHelpers.h"
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
  static constexpr std::string_view cutNames[4] = {"DiffCuts2pi", "DiffCuts4pi", "DiffCuts2K", "DiffCuts4K"};
  MutableConfigurable<DGCutparHolder> diffCuts2pi{cutNames[0].data(), {}, "Diffractive 2pi events cuts"};
  MutableConfigurable<DGCutparHolder> diffCuts4pi{cutNames[1].data(), {}, "Diffractive 4pi events cuts"};
  MutableConfigurable<DGCutparHolder> diffCuts2K{cutNames[2].data(), {}, "Diffractive 2K events cuts"};
  MutableConfigurable<DGCutparHolder> diffCuts4K{cutNames[3].data(), {}, "Diffractive 4K events cuts"};

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
  static constexpr std::string_view histNames[4] = {"aftercut2pi", "aftercut4pi", "aftercut2K", "aftercut4K"};
  HistogramRegistry registry{
    "registry",
    {
      {histNames[0].data(), "#aftercut2pi", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
      {histNames[1].data(), "#aftercut4pi", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
      {histNames[2].data(), "#aftercut2K", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
      {histNames[3].data(), "#aftercut4K", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
    }};

  void init(InitContext&)
  {
  }

  // some general Collisions and Tracks filter
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

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

    // nominal BC
    auto bc = collision.bc_as<BCs>();

    // loop over 4 cases
    bool ccs[4]{false};
    for (int ii = 0; ii < 4; ii++) {

      // different cases
      switch (ii) {
        case 0:
          diffCuts = (DGCutparHolder)diffCuts2pi;
          registry.fill(HIST(histNames[0]), 0.);
          break;
        case 1:
          diffCuts = (DGCutparHolder)diffCuts4pi;
          registry.fill(HIST(histNames[1]), 0.);
          break;
        case 2:
          diffCuts = (DGCutparHolder)diffCuts2K;
          registry.fill(HIST(histNames[2]), 0.);
          break;
        case 3:
          diffCuts = (DGCutparHolder)diffCuts4K;
          registry.fill(HIST(histNames[3]), 0.);
          break;
        default:
          continue;
      }

      // obtain slice of compatible BCs
      auto bcRange = compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
      LOGF(debug, "  Number of compatible BCs in +- %i / %i dtcoll: %i", diffCuts.NDtcoll(), diffCuts.minNBCs(), bcRange.size());

      // apply DG selection
      auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bc, bcRange, tracks, fwdtracks);

      // save decision
      if (isDGEvent == 0) {
        LOGF(debug, "This collision is a DG candidate!");
      }
      ccs[ii] = (isDGEvent == 0);

      // update after cut histogram
      // different cases
      switch (ii) {
        case 0:
          registry.fill(HIST(histNames[0]), isDGEvent + 1);
          break;
        case 1:
          registry.fill(HIST(histNames[1]), isDGEvent + 1);
          break;
        case 2:
          registry.fill(HIST(histNames[2]), isDGEvent + 1);
          break;
        case 3:
          registry.fill(HIST(histNames[3]), isDGEvent + 1);
          break;
        default:
          continue;
      }
    }

    // update filterTable
    filterTable(ccs[0], ccs[1], ccs[2], ccs[3]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGFilterRun3>(cfgc, TaskName{"DGfilterRun3"}),
  };
}
