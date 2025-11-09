// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \file phitutorial.cxx
/// \brief Phi meson analysis tutorial
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

// IMPORTANT INCLUDES
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/ASoAHelpers.h>
#include <Framework/runDataProcessing.h>

// ROOT Includes (optional)
#include <TLorentzVector.h>

// C++ includes
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// MAIN STRUCT
struct phitutorial_step0 {

  //*************************************//
  // SLICECACHE AND REGISTRY DEFS
  //*************************************//
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //*************************************//
  // INIT FUNCTION AND HISTOGRAM BOOKING
  //*************************************//
  void init(o2::framework::InitContext&)
  {
    const AxisSpec ptAxis = {200, 0, 20.0};
    const AxisSpec MinvAxis = {200, 0.85, 1.25};

    histos.add("Nch_pT", "Nch_pT", kTH1F, {ptAxis});

  }; // end of init

  //*************************************//
  // TIME TO BUILD TRACK AND EVENT CANDIDATES
  //*************************************//
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::CentFT0Ms>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  double massKa = o2::constants::physics::MassKPlus;

  //***************************************//
  // PREAMBLE COMPLETE, NOW WE DO HELPER FCNS
  //**************************************//
  template <typename EventType>
  bool eventSelection(const EventType event)
  {
    if (!event.sel8()) // This is required to extract good events
      return false;

    return true;
  };
  //********************************************//
  // Space for more helper functions!

  //********************************************//
  // HELPER FCNS COMPLETE, NOW WE DO PROCESS FCNS
  //********************************************//

  // SAME EVENT
  int nEvents = 0;
  void processDataSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    nEvents++;
    if ((nEvents + 1) % 10000 == 0) {
      std::cout << "Processed Data Events: " << nEvents << std::endl;
    }

    if (!eventSelection(collision))
      return;

    // Now, time to start coding the task!
    // Keep in mind that:
    // M_inv = sqrt( (E1+E2)^2 - |P1 + P2|^2 )
    // Where you use the energies and momenta of the individual Kaons.

    // You should fill: histos.fill(HIST("Minv"), M_inv), calculated as above.

    // Usefull tips:
    //  E = sqrt(p^2 + m^2). The Kaon mass is found above in the constant massKa
    //  pz = pT*sinh(eta)
    //  track.pt()
    //  track.eta();
    //  std::sinh(x)

    // For more concise techinques, check out:
    // ROOT::Math::PxPyPzMVector //Check google
    // combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy....   //check ALICE O2 documentation

    for (const auto& track : tracks) {
      histos.fill(HIST("Nch_pT"), track.pt());
      //..
      //..
      //..
    }

  } // proccessSameEvent
  PROCESS_SWITCH(phitutorial_step0, processDataSameEvent, "process Data Same Event", false);

  //***************************************//
  // TASK COMPLETE!
  //**************************************//
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phitutorial_step0>(cfgc)};
};
