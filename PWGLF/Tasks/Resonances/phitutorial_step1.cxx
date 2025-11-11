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
struct phitutorial_step1 {

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
    histos.add("Nch_USS_Minv", "Nch_USS_Minv", kTH1F, {MinvAxis});

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

    // Now, we want to add some kind of filter for our tracks!
    // Tracks we want to accept:
    // track.isGlobalTrack()  <.... this menas that it is a good track
    // track.pt() >0.15       <.... we want to remove really low momentum tracks
    // -0.8<track.eta()<0.8   <.... the TPC efficiency drops really hard at the edges

    // We can either create a template function (see eventSelection above for reference) that accepts tracks instead of collisions.
    // We can also directly implement the track cuts on the track iterators (track, trk1, trk2) themselves.
    // We can also set a filter around the subscribed table. Then before the process function, you would have to modify:
    // Filter pT Filter = track::pt > 0.15f;
    // then you have to modify the subscription
    // soa::Filtered<TrackCandidates> const& tracks

    for (const auto& track : tracks) {
      histos.fill(HIST("Nch_pT"), track.pt());
    }

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

      ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);

      lResonance = lDecayDaughter1 + lDecayDaughter2;
      double conjugate = trk1.sign() * trk2.sign();
      if (conjugate < 0) {
        histos.fill(HIST("Nch_USS_Minv"), lResonance.M());
      }
    } // Invariant mass combinations

  } // proccessSameEvent
  PROCESS_SWITCH(phitutorial_step1, processDataSameEvent, "process Data Same Event", false);

  //***************************************//
  // TASK COMPLETE!
  //**************************************//
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phitutorial_step1>(cfgc)};
};
