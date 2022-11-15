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

#include "Framework/Logger.h"
#include "TLorentzVector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/DGSelector.h"

// -----------------------------------------------------------------------------
DGSelector::DGSelector()
{
  fPDG = TDatabasePDG::Instance();
}

// -----------------------------------------------------------------------------
DGSelector::~DGSelector()
{
  delete fPDG;
}

// -----------------------------------------------------------------------------
// Function to check if collisions passes filter
template <typename CC, typename BCs, typename TCs, typename FWs>
int DGSelector::IsSelected(DGCutparHolder diffCuts, CC const& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
{
  LOGF(debug, "Collision %f", collision.collisionTime());
  LOGF(debug, "Number of close BCs: %i", bcRange.size());

  // check that there are no FIT signals in any of the compatible BCs
  // Double Gap (DG) condition
  auto lims = diffCuts.FITAmpLimits();

  for (auto const& bc : bcRange) {
    LOGF(debug, "Amplitudes FV0A %f FT0 %f / %f FDD %i / %i",
         bc.has_foundFV0() ? FV0AmplitudeA(bc.foundFV0()) : -1.,
         bc.has_foundFT0() ? FT0AmplitudeA(bc.foundFT0()) : -1.,
         bc.has_foundFT0() ? FT0AmplitudeC(bc.foundFT0()) : -1.,
         bc.has_foundFDD() ? FDDAmplitudeA(bc.foundFDD()) : -1,
         bc.has_foundFDD() ? FDDAmplitudeC(bc.foundFDD()) : -1);
    LOGF(debug, "  clean FV0A %i FT0 %i FDD %i", cleanFV0(bc, lims[0]), cleanFT0(bc, lims[1], lims[2]), cleanFDD(bc, lims[3], lims[4]));

    if (!udhelpers::cleanFIT(bc, diffCuts.FITAmpLimits())) {
      return 1;
    }
  }

  // no activity in muon arm
  LOGF(debug, "FwdTracks %i", fwdtracks.size());
  for (auto& fwdtrack : fwdtracks) {
    LOGF(debug, "  %i / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p());
  }
  if (fwdtracks.size() > 0) {
    return 2;
  }

  // no global tracks which are not vtx tracks
  // no vtx tracks which are not global tracks
  auto rgtrwTOF = 0.;
  for (auto& track : tracks) {
    if (track.isGlobalTrack() && !track.isPVContributor()) {
      return 3;
    }
    if (diffCuts.globalTracksOnly() && !track.isGlobalTrack() && track.isPVContributor()) {
      return 4;
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
    return 5;
  }

  // number of vertex tracks
  if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
    return 6;
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
      if (!hasGoodPID(diffCuts, track)) {
        return 7;
      }

      // pt
      lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
      if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
        return 8;
      }

      // eta
      if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
        return 9;
      }
      netCharge += track.sign();
      ivm += lvtmp;
    }
  }

  // net charge
  auto netChargeValues = diffCuts.netCharges();
  if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
    return 10;
  }
  // invariant mass
  if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
    return 11;
  }

  // if we arrive here then the event is good!
  return 0;
};

// -----------------------------------------------------------------------------
/*
template <typename BCs, typename TCs, typename FWs>
int DGSelector::IsSelected(DGCutparHolder diffCuts, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
{
  // check that there are no FIT signals in bcRange
  // Double Gap (DG) condition
  for (auto const& bc : bcRange) {
    if (!udhelpers::cleanFIT(bc, diffCuts.FITAmpLimits())) {
      return 1;
    }
  }

  // no activity in muon arm
  LOGF(debug, "FwdTracks %i", fwdtracks.size());
  for (auto& fwdtrack : fwdtracks) {
    LOGF(debug, "  %i / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p());
  }
  if (fwdtracks.size() > 0) {
    return 2;
  }

  // number of tracks
  if ((int)tracks.size() < diffCuts.minNTracks() || (int)tracks.size() > diffCuts.maxNTracks()) {
    return 6;
  }

  // PID, pt, and eta of tracks, invariant mass, and net charge
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
    if (!hasGoodPID(diffCuts, track)) {
      return 7;
    }

    // pt
    lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
    if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
      return 8;
    }

    // eta
    if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
      return 9;
    }
    netCharge += track.sign();
    ivm += lvtmp;
  }

  // net charge
  auto netChargeValues = diffCuts.netCharges();
  if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
    return 10;
  }

  // invariant mass
  if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
    return 11;
  }

  // if we arrive here then the event is good!
  return 0;
};
*/
// -----------------------------------------------------------------------------
