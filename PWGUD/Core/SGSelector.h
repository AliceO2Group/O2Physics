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

#ifndef PWGUD_CORE_SGSELECTOR_H_
#define PWGUD_CORE_SGSELECTOR_H_

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGCutParHolder.h"

// Selector for Single Gap events
class SGSelector
{
 public:
  // constructor/destructor
  SGSelector() { fPDG = TDatabasePDG::Instance(); }
  ~SGSelector() { delete fPDG; }

  template <typename CC, typename BCs, typename TCs, typename FWs>
  int Print(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(info, "Size of array %i", collision.size());
    return 1;
  }

  // Function to check if collisions pass SG filter
  template <typename CC, typename BCs, typename TCs, typename FWs>
  int IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(debug, "Collision %f", collision.collisionTime());
    LOGF(debug, "Number of close BCs: %i", bcRange.size());
    // check that there are no FIT signals in any of the compatible BCs
    // Single Gap (SG) condition
    bool gA = true;
    bool gC = true;
    for (auto const& bc : bcRange) {
      if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        gA = false;
      }
      if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        gC = false;
      }
    }
    if (!gA && !gC)
      return 3;
    // forward tracks
    LOGF(debug, "FwdTracks %i", fwdtracks.size());
    if (!diffCuts.withFwdTracks()) {
      for (auto& fwdtrack : fwdtracks) {
        LOGF(debug, "  %i / %f / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p(), fwdtrack.trackTimeRes());
        // only consider tracks with MID (good timing)
        if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) {
          return 4;
        }
      }
    }

    // no global tracks which are not vtx tracks
    // no vtx tracks which are not global tracks
    // no PV tracks with ITS only
    auto rgtrwTOF = 0.; // fraction of PV tracks with TOF hit
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        return 5;
      }
      if (diffCuts.globalTracksOnly() && track.isPVContributor() && !track.isGlobalTrack()) {
        return 6;
      }
      if (!diffCuts.ITSOnlyTracks() && track.isPVContributor() && !track.hasTPC()) {
        return 7;
      }

      // update fraction of PV tracks with TOF hit
      if (track.isPVContributor() && track.hasTOF()) {
        rgtrwTOF += 1.;
      }
    }
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    if (rgtrwTOF < diffCuts.minRgtrwTOF()) {
      return 8;
    }

    // number of vertex tracks
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      return 9;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // consider only vertex tracks

    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      if (track.isPVContributor()) {

        // PID
        // if (!udhelpers::hasGoodPID(diffCuts, track)) {
        //   return 8;
        // }

        // pt
        lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
        if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
          return 10;
        }

        // eta
        if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
          return 11;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }

    // net charge
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      return 12;
    }
    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 13;
    }
    // Single Gap (SG) condition
    if (gA && gC) {
      return 2;
    } else if (gA) {
      return 0;
    } else if (gC) {
      return 1;
    }
  };

  // Function to check if BC passes SG filter (without associated collision)
  template <typename BCs, typename TCs, typename FWs>
  int IsSelected(SGCutParHolder diffCuts, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    // check that there are no FIT signals in bcRange
    // Single Gap (SG) condition
    bool gA = true;
    bool gC = true;
    for (auto const& bc : bcRange) {
      if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        gA = false;
      }
      if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        gC = false;
      }
    }
    if (!gA && !gC)
      return 3;
    // no activity in muon arm
    if (!diffCuts.withFwdTracks()) {
      for (auto& fwdtrack : fwdtracks) {
        LOGF(debug, "  %i / %f / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p(), fwdtrack.trackTimeRes());
        // only consider tracks with MID (good timing)
        if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) {
          return 4;
        }
      }
    }

    // number of tracks
    if (static_cast<int>(tracks.size()) < diffCuts.minNTracks() || static_cast<int>(tracks.size()) > diffCuts.maxNTracks()) {
      return 5;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // consider only vertex tracks

    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {

      // PID
      // if (!udhelpers::hasGoodPID(diffCuts, track)) {
      //   return 4;
      // }

      // pt
      lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
      if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
        return 6;
      }

      // eta
      if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
        return 7;
      }
      netCharge += track.sign();
      ivm += lvtmp;
    }

    // net charge
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      return 8;
    }
    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 9;
    }

    // check that there are no FIT signals in any of the compatible BCs
    // Single Gap (SG) condition
    if (gA && gC) {
      return 2;
    } else if (gA) {
      return 0;
    } else if (gC) {
      return 1;
    }
  };

 private:
  TDatabasePDG* fPDG;
};

#endif // PWGUD_CORE_SGSELECTOR_H_
