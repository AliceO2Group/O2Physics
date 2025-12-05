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
/// \file   tofSkimsTableCreator.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  19/11/2022
/// \brief  Task to defined the skimmed data format for the TOF skims
///

#include "tofSkimsTableCreator.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

struct tofSkimsTableCreator {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TOFEvTime, aod::EvTimeTOFOnly, aod::TOFSignal, aod::pidEvTimeFlags,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>;

  // Tables to be produced
  Produces<o2::aod::SkimmedTOFColl> tableColRow;
  Produces<o2::aod::SkimmedTOF> tableRow;

  // Configurables
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<bool> keepTpcOnly{"keepTpcOnly", false, "Flag to keep the TPC only tracks as well"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 0.1, "Fractions of events to keep"};

  unsigned int randomSeed = 0;
  void init(o2::framework::InitContext&)
  {
    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    switch (applyEvSel.value) {
      case 0:
      case 1:
      case 2:
        break;
      default:
        LOG(fatal) << "Invalid event selection flag: " << applyEvSel.value;
        break;
    }
    switch (trackSelection.value) {
      case 0:
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
        break;
      default:
        LOG(fatal) << "Invalid track selection flag: " << trackSelection.value;
        break;
    }
  }

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (o2::aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (o2::aod::evsel::sel8 == true));
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireInAcceptanceTracksInFilter());

  void process(soa::Filtered<Coll>::iterator const& collision,
               soa::Filtered<Trks> const& tracks)
  {
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return;
    }
    tableRow.reserve(tracks.size());
    float evTimeT0AC = 0.f;
    float evTimeT0ACErr = 0.f;
    if (collision.t0ACValid()) {
      evTimeT0AC = collision.t0AC() * 1000.f;
      evTimeT0ACErr = collision.t0resolution() * 1000.f;
    }
    float evTimeT0A = 0.f;
    if (collision.t0ACorrectedValid()) {
      evTimeT0A = collision.t0ACorrected() * 1000.f;
    }
    float evTimeT0C = 0.f;
    if (collision.t0CCorrectedValid()) {
      evTimeT0C = collision.t0CCorrected() * 1000.f;
    }
    float evTimeTOF = 0.f;
    float evTimeTOFErr = 0.f;
    int evTimeTOFMult = 0.f;
    uint8_t tofFlags = 0;
    for (auto const& trk : tracks) {
      if (!trk.hasTOF()) {
        continue;
      }
      evTimeTOF = trk.evTimeTOF();
      evTimeTOFErr = trk.evTimeTOFErr();
      evTimeTOFMult = trk.evTimeTOFMult();
      tofFlags = trk.tofFlags();
      break;
    }

    tableColRow(evTimeTOF,
                evTimeTOFErr,
                evTimeTOFMult,
                evTimeT0A,
                evTimeT0C,
                evTimeT0AC,
                evTimeT0ACErr,
                collision.collisionTime(),
                collision.collisionTimeRes(),
                tofFlags);

    int8_t lastTRDLayer = -1;
    for (auto const& trk : tracks) {
      if (!keepTpcOnly.value && !trk.hasTOF()) {
        continue;
      }

      lastTRDLayer = -1;
      if (trk.hasTRD()) {
        for (int8_t l = 7; l >= 0; l--) {
          if (trk.trdPattern() & (1 << l)) {
            lastTRDLayer = l;
            break;
          }
        }
      }
      tableRow(tableColRow.lastIndex(),
               trk.p(),
               trk.pt() * trk.sign(),
               trk.eta(),
               trk.phi(),
               trk.pidForTracking(),
               trk.tofExpMom(),
               trk.length(),
               trk.tofChi2(),
               trk.tpcSignal(),
               trk.tofSignal(),
               trk.evTimeTOF(),
               trk.evTimeTOFErr(),
               evTimeT0AC,
               evTimeT0ACErr,
               trk.tofFlags(),
               trk.tpcInnerParam(),
               trk.tpcNClsFindable(),
               trk.tpcNClsFindableMinusFound(),
               trk.tpcNClsFindableMinusCrossedRows(),
               trk.tpcNClsShared(),
               lastTRDLayer);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofSkimsTableCreator>(cfgc)}; }
