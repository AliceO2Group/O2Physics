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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
/// O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include "tofSkimsTableCreator.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;

struct tofSimsTableCreator {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::TOFEvTime, aod::EvTimeTOFOnly, aod::TOFSignal, aod::pidEvTimeFlags,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::FT0sCorrected>;

  // Tables to be produced
  Produces<o2::aod::SkimmedTOF> tableRow;

  // Configurables
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> applyTrkSel{"applyTrkSel", 1, "Flag to apply track selection: 0 -> no track selection, 1 -> track selection"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 0.1, "Fractions of events to keep"};

  void init(o2::framework::InitContext& initContext) {}

  void process(Coll::iterator const& collision,
               Trks const& tracks)
  {
    if (fractionOfEvents < 1.f && (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return;
    }

    switch (applyEvSel.value) {
      case 0:
        break;
      case 1:
        if (!collision.sel7()) {
          return;
        }
        break;
      case 2:
        if (!collision.sel8()) {
          return;
        }
        break;
      default:
        LOG(fatal) << "Invalid event selection flag: " << applyEvSel.value;
        break;
    }
    tableRow.reserve(tracks.size());
    float evTimeT0AC = 0.f;
    float evTimeT0ACErr = 0.f;
    if (collision.t0ACValid()) {
      evTimeT0AC = collision.t0AC() * 1000.f;
      evTimeT0ACErr = collision.t0resolution() * 1000.f;
    }

    for (auto const& trk : tracks) {
      tableRow(trk.collisionId(),
               trk.p(),
               trk.pt(),
               trk.eta(),
               trk.phi(),
               trk.pidForTracking(),
               trk.tofExpMom(),
               trk.length(),
               trk.tofChi2(),
               trk.tofSignal(),
               trk.evTimeTOF(),
               trk.evTimeTOFErr(),
               evTimeT0AC,
               evTimeT0ACErr,
               trk.tofFlags());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofSimsTableCreator>(cfgc)};
}
