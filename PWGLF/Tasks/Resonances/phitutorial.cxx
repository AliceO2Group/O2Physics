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
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
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
struct phitutorial {

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
    histos.add("Nch_LSS_Minv", "Nch_LSS_Minv", kTH1F, {MinvAxis});

    histos.add("Nch_ME_Minv", "Nch_ME_Minv", kTH1F, {MinvAxis});

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
    if (!event.sel8())
      return false;
    if (std::abs(event.posZ()) > 10)
      return false;
    if (!event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (!event.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (!event.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
      return false;

    return true;
  };
  //********************************************//
  template <typename TracksType>
  bool trackSelection(const TracksType track)
  {
    if (!track.isGlobalTrack())
      return false;
    if (track.pt() < 0.15)
      return false;
    if (std::abs(track.eta()) > 1.0)
      return false;

    return true;
  };
  //********************************************//
  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    // TPC
    if (std::abs(candidate.tpcNSigmaKa()) < 3)
      tpcPIDPassed = true;
    // TOF
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < 3) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

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

    for (const auto& track : tracks) {
      histos.fill(HIST("Nch_pT"), track.pt());
    }

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      if (!trackSelection(trk1) || !trackSelection(trk2)) {
        continue;
      }
      if (!trackPIDKaon(trk1) || !trackPIDKaon(trk2)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);

      lResonance = lDecayDaughter1 + lDecayDaughter2;
      double conjugate = trk1.sign() * trk2.sign();
      if (conjugate < 0) {
        histos.fill(HIST("Nch_USS_Minv"), lResonance.M());
      } else {
        histos.fill(HIST("Nch_LSS_Minv"), lResonance.M());
      }
    } // Invariant mass combinations

  } // proccessSameEvent
  PROCESS_SWITCH(phitutorial, processDataSameEvent, "process Data Same Event", false);

  //**************************************************************************************************************************//

  // MIXED EVENT

  //*********************************************************//
  // DEFINITION OF SLICE CACHE, BINNING AND MIXING STRUCTURE
  //*********************************************************//
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  std::vector<double> zBins{10, -10, 10};
  std::vector<double> multBins{VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType binning{{zBins, multBins}, true};
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binning, 5, -1, &cache};

  void processDataMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    for (const auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!eventSelection(c1) || !eventSelection(c2))
        continue;

      for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (track1.sign() * track2.sign() > 0)
          continue;

        if (!trackSelection(track1) || !trackSelection(track2)) {
          continue;
        }
        if (!trackPIDKaon(track1) || !trackPIDKaon(track2)) {
          continue;
        }

        ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, mother;
        lDecayDaughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        lDecayDaughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);

        mother = lDecayDaughter1 + lDecayDaughter2;

        histos.fill(HIST("Nch_ME_Minv"), mother.M());
      }
    }
  } // processMixedEvent
  PROCESS_SWITCH(phitutorial, processDataMixedEvent, "process Data Mixed Event", false);
};

//***************************************//
// TASK COMPLETE!
//**************************************//

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phitutorial>(cfgc)};
};
