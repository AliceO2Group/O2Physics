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

#ifndef O2_ANALYSIS_DIFFRACTION_SELECTOR_H_
#define O2_ANALYSIS_DIFFRACTION_SELECTOR_H_

#include "Framework/DataTypes.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename TC>
bool hasGoodPID(cutHolder diffCuts, TC track);

// add here Selectors for different types of diffractive events
// Selector for Double Gap events
struct DGSelector {
 public:
  DGSelector() = default;

  // Function to check if collisions passes filter
  template <typename CC, typename BC, typename BCs, typename TCs, typename FWs>
  bool IsSelected(cutHolder diffCuts, CC const& collision, BC& bc, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(debug, "Collision %f BC %i", collision.collisionTime(), bc.globalBC());
    LOGF(debug, "Number of close BCs: %i", bcRange.size());
    LOGF(debug, "Number of tracks: %i", tracks.size());

    // Number of tracks
    // All tracks in vertex
    if (collision.numContrib() != tracks.size()) {
      return false;
    }
    if (tracks.size() < diffCuts.minNTracks() || tracks.size() > diffCuts.maxNTracks()) {
      return false;
    }

    // all tracks must be global tracks
    // and some need to have a TOF hit
    // chaeck also net charge
    auto goodTracks = true;
    int nTracksWithTOFHit = 0;
    auto netCharge = 0;
    LOGF(debug, "Tracks");
    for (auto& track : tracks) {
      goodTracks &= (track.isGlobalTrack() > 0);

      // update netCharge
      netCharge += track.sign();

      // TOF hit with good match quality
      if (track.hasTOF() && track.tofChi2() < diffCuts.maxTOFChi2()) {
        nTracksWithTOFHit++;
      }
      LOGF(debug, "   global %i, TOF [%i] signal %f / chi2: %f", track.isGlobalTrack(), (track.hasTOF() ? 1 : 0), track.tofSignal(), track.tofChi2());
    }
    LOGF(debug, "  Tracks net charge %i, with TOF hit %i / %i", netCharge, nTracksWithTOFHit, diffCuts.minNTracksWithTOFHit());
    goodTracks &= (nTracksWithTOFHit >= diffCuts.minNTracksWithTOFHit());
    goodTracks &= (netCharge >= diffCuts.minNetCharge() && netCharge <= diffCuts.maxNetCharge());
    if (!goodTracks) {
      return false;
    }

    // only tracks with good PID
    auto goodPID = true;
    for (auto& track : tracks) {
      goodPID &= hasGoodPID(diffCuts, track);
    }
    if (!goodPID) {
      return false;
    }

    // check no activity in muon arm
    LOGF(debug, "Muons %i", fwdtracks.size());
    for (auto& muon : fwdtracks) {
      LOGF(debug, "  %i / %f / %f / %f", muon.trackType(), muon.eta(), muon.pt(), muon.p());
    }
    if (fwdtracks.size() > 0) {
      return false;
    }

    // ATTENTION: currently all events are selected
    return true;
  };
};

// function to check if track provides good PID information
// Checks the nSigma for any particle assumption to be within limits.
template <typename TC>
bool hasGoodPID(cutHolder diffCuts, TC track)
{
  // El, Mu, Pi, Ka, and Pr are considered
  // at least one nSigma must be within set limits
  LOGF(debug, "TPC PID %f / %f / %f / %f / %f",
       track.tpcNSigmaEl(),
       track.tpcNSigmaMu(),
       track.tpcNSigmaPi(),
       track.tpcNSigmaKa(),
       track.tpcNSigmaPr());
  if (TMath::Abs(track.tpcNSigmaEl()) < diffCuts.maxnSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaMu()) < diffCuts.maxnSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaPi()) < diffCuts.maxnSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaKa()) < diffCuts.maxnSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaPr()) < diffCuts.maxnSigmaTPC()) {
    return true;
  }

  if (track.hasTOF()) {
    LOGF(debug, "TOF PID %f / %f / %f / %f / %f",
         track.tofNSigmaEl(),
         track.tofNSigmaMu(),
         track.tofNSigmaPi(),
         track.tofNSigmaKa(),
         track.tofNSigmaPr());
    if (TMath::Abs(track.tofNSigmaEl()) < diffCuts.maxnSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaMu()) < diffCuts.maxnSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaPi()) < diffCuts.maxnSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaKa()) < diffCuts.maxnSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaPr()) < diffCuts.maxnSigmaTOF()) {
      return true;
    }
  }
  return false;
}

#endif // O2_ANALYSIS_DIFFRACTION_SELECTOR_H_
