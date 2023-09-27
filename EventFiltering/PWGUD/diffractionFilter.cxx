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
#include "Framework/StaticFor.h"
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

  // Configurables
  int64_t collCounter;
  int64_t verbose;
  Configurable<int64_t> Verbose{"Verbose", {0}, "With extra log output"};
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> diffCutsHolder{"DiffCuts", {}, "Diffractive events cuts"};

  // DG selector
  DGSelector dgSelector;

  // histogram stat/aftercuts with cut statistics
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
  // initinitialize HistogramRegistry
  static constexpr std::string_view hcFITs[5] = {"FIT/cleanFITFV0A", "FIT/cleanFITFT0A", "FIT/cleanFITFT0C", "FIT/cleanFITFDDA", "FIT/cleanFITFDDC"};
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // reset collision counter
    collCounter = 0;
    verbose = (int64_t)Verbose;
    if (verbose <= 0) {
      verbose = 1000000000000;
    }

    // cut holder
    diffCuts = (DGCutparHolder)diffCutsHolder;

    // create histograms
    // stat
    registry.add("stat/aftercuts", "Cut efficiencies", {HistType::kTH1F, {{14, -0.5, 13.5}}});

    // FIT
    registry.add("FIT/cleanFIT", "Rejection by FIT veto versus tested BC range", {HistType::kTH2F, {{21, -0.5, 20.5}, {2, -0.5, 1.5}}});
    for (auto n{0}; n < 5; n++) {
      registry.add(hcFITs[n].data(), hcFITs[n].data(), {HistType::kTH2F, {{21, -0.5, 20.5}, {2, -0.5, 1.5}}});
    }
    registry.add("FIT/FV0Atime", "Time of FV0A [ns]", {HistType::kTH1F, {{6001, -20.005, 40.005}}});
    registry.add("FIT/FT0Atime", "Time of FT0A [ns]", {HistType::kTH1F, {{6001, -20.005, 40.005}}});
    registry.add("FIT/FT0Ctime", "Time of FT0C [ns]", {HistType::kTH1F, {{6001, -20.005, 40.005}}});
    registry.add("FIT/FDDAtime", "Time of FDDA [ns]", {HistType::kTH1F, {{6001, -20.005, 40.005}}});
    registry.add("FIT/FDDCtime", "Time of FDDC [ns]", {HistType::kTH1F, {{6001, -20.005, 40.005}}});

    // collision
    registry.add("collisions/tracksAll", "Number of tracks per collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/PVTracksAll", "Number of PV tracks per collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/globalTracksAll", "Number of global tracks per collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/netChargeAll", "Net charge of PV tracks versus number of PV tracks of all collisions", {HistType::kTH2F, {{51, -0.5, 50.5}, {41, -20.5, 20.5}}});
    registry.add("collisions/dtcvsrPVtrwTOFAll", "Collision time resolution as function of the fraction of PV tracks with TOF hit of all collisions", {HistType::kTH2F, {{1001, -0.05, 100.05}, {101, -0.005, 1.005}}});
    registry.add("collisions/rPVtrwTOFAll", "Fraction of PV tracks with TOF hit versus number of PV tracks of all collisions", {HistType::kTH2F, {{100, 0.5, 100.5}, {101, -0.005, 1.005}}});
    registry.add("collisions/forwardTracksAll", "Number of forward tracks per collision", {HistType::kTH2F, {{5, -0.5, 4.5}, {101, -0.5, 100.5}}});

    registry.add("collisions/tracksDG", "Number of tracks per DG collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/PVTracksDG", "Number of PV tracks per DG collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/globalTracksDG", "Number of global tracks per DG collision", {HistType::kTH1F, {{300, 0.5, 300.5}}});
    registry.add("collisions/netChargeDG", "Net charge of PV tracks versus number of PV tracks of DG collisions", {HistType::kTH2F, {{51, -0.5, 50.5}, {41, -20.5, 20.5}}});
    registry.add("collisions/dtcvsrPVtrwTOFDG", "Collision time resolution as function of the fraction of PV tracks with TOF hit of DG candidates", {HistType::kTH2F, {{1001, -0.05, 100.05}, {101, -0.005, 1.005}}});
    registry.add("collisions/rPVtrwTOFDG", "Fraction of PV tracks with TOF hit versus number of PV tracks of DG collisions", {HistType::kTH2F, {{100, 0.5, 100.5}, {101, -0.005, 1.005}}});
    registry.add("collisions/forwardTracksDG", "Number of forward tracks per DG collision", {HistType::kTH2F, {{5, -0.5, 4.5}, {101, -0.5, 100.5}}});

    // tracks
    registry.add("tracks/etavsptAll", "eta versus pt of PV tracks of all collisions", {HistType::kTH2F, {{401, -2.005, 2.005}, {1000, 0., 10.}}});
    registry.add("tracks/etavsptDGwT", "eta versus pt of PV tracks with TOF hit in DG collisions", {HistType::kTH2F, {{401, -2.005, 2.005}, {1000, 0., 10.}}});
    registry.add("tracks/etavsptDGnT", "eta versus pt of PV tracks without TOF hit in DG collisions", {HistType::kTH2F, {{401, -2.005, 2.005}, {1000, 0., 10.}}});

    // forwardTracks
    registry.add("forwardTracks/timeResolution", "Time resolution of forward tracks as function of track type [ns]", {HistType::kTH2F, {{5, -0.5, 4.5}, {10001, -0.005, 100.005}}});
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
    registry.fill(HIST("stat/aftercuts"), 0.);

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());
    LOGF(debug, "  Number of compatible BCs in +- %i / %i dtcoll: %i", diffCuts.NDtcoll(), diffCuts.minNBCs(), bcRange.size());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bcRange, tracks, fwdtracks);

    // update after cut histogram
    registry.fill(HIST("stat/aftercuts"), isDGEvent + 1);

    // update filterTable
    ccs = (isDGEvent == 0);
    filterTable(ccs);

    // log output to check consistency of selections on original and skimmed data
    auto bc2 = collision.foundBC_as<BCs>();
    auto rgtrwTOF = udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib());
    if ((collCounter % verbose) == 0) {
      auto bc1 = collision.bc_as<BCs>();
      LOGF(info, "BCId %d/%d isDG %d PV tracks %d rgtrwTOF %f", bc1.globalBC(), bc2.globalBC(), isDGEvent + 1, collision.numContrib(), rgtrwTOF);
    }

    // update QC histograms
    // FIT
    if (bc2.has_foundFV0()) {
      registry.fill(HIST("FIT/FV0Atime"), bc2.foundFV0().time());
    }
    if (bc2.has_foundFT0()) {
      registry.fill(HIST("FIT/FT0Atime"), bc2.foundFT0().timeA());
      registry.fill(HIST("FIT/FT0Ctime"), bc2.foundFT0().timeC());
    }
    if (bc2.has_foundFDD()) {
      registry.fill(HIST("FIT/FDDAtime"), bc2.foundFDD().timeA());
      registry.fill(HIST("FIT/FDDCtime"), bc2.foundFDD().timeC());
    }

    auto FITlims = std::vector<float>(5, 1000000.);
    bool isDGcandidate = true;
    for (int nMinBC = 0; nMinBC <= 20; nMinBC++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, 0, bcs, nMinBC);
      isDGcandidate = true;
      for (auto const& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());
      }
      registry.get<TH2>(HIST("FIT/cleanFIT"))->Fill(nMinBC, isDGcandidate * 1.);

      // loop over single detectors
      static_for<0, 4>([&](auto n) {
        FITlims[n] = 0.;
        isDGcandidate = true;
        for (auto const& bc : bcSlice) {
          isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), FITlims);
        }
        constexpr int index = n.value;
        registry.fill(HIST(hcFITs[index]), nMinBC, isDGcandidate * 1.);
        FITlims[n] = 1000000.;
      });
    }

    // collisions
    registry.fill(HIST("collisions/tracksAll"), tracks.size());
    registry.fill(HIST("collisions/PVTracksAll"), collision.numContrib());
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    registry.get<TH1>(HIST("collisions/globalTracksAll"))->Fill(goodTracks.size());
    auto netCharge = udhelpers::netCharge<true>(tracks);
    registry.fill(HIST("collisions/netChargeAll"), collision.numContrib(), netCharge);
    registry.fill(HIST("collisions/dtcvsrPVtrwTOFAll"), collision.collisionTimeRes(), rgtrwTOF);
    registry.fill(HIST("collisions/rPVtrwTOFAll"), collision.numContrib(), rgtrwTOF);
    if (ccs) {
      registry.fill(HIST("collisions/tracksDG"), tracks.size());
      registry.fill(HIST("collisions/PVTracksDG"), collision.numContrib());
      registry.get<TH1>(HIST("collisions/globalTracksDG"))->Fill(goodTracks.size());
      registry.fill(HIST("collisions/netChargeDG"), collision.numContrib(), netCharge);
      registry.fill(HIST("collisions/dtcvsrPVtrwTOFDG"), collision.collisionTimeRes(), rgtrwTOF);
      registry.fill(HIST("collisions/rPVtrwTOFDG"), collision.numContrib(), rgtrwTOF);
    }

    // tracks
    for (auto const& track : tracks) {
      if (track.isPVContributor()) {
        registry.fill(HIST("tracks/etavsptAll"), track.eta(), track.pt());
        if (ccs) {
          if (track.hasTOF()) {
            registry.fill(HIST("tracks/etavsptDGwT"), track.eta(), track.pt());
          } else {
            registry.fill(HIST("tracks/etavsptDGnT"), track.eta(), track.pt());
          }
        }
      }
    }

    // forward tracks
    int nforwardTracks[5] = {0};
    for (auto track : fwdtracks) {
      nforwardTracks[track.trackType()]++;
      registry.fill(HIST("forwardTracks/timeResolution"), track.trackType(), track.trackTimeRes());
    }
    for (auto ii = 0; ii < 5; ii++) {
      registry.fill(HIST("collisions/forwardTracksAll"), ii, nforwardTracks[ii]);
      if (ccs) {
        registry.fill(HIST("collisions/forwardTracksDG"), ii, nforwardTracks[ii]);
      }
    }

    // update collision counter
    collCounter++;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGFilterRun3>(cfgc, TaskName{"DGfilterRun3"}),
  };
}
