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

class SGSelector
{
 public:
  SGSelector() : fPDG(TDatabasePDG::Instance()) {}

  template <typename CC, typename BCs, typename TCs, typename FWs>
  int Print(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(info, "Size of array %i", collision.size());
    return 1;
  }

  template <typename CC, typename BCs>
  int IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange)
  {
    LOGF(debug, "Collision %f", collision.collisionTime());
    LOGF(debug, "Number of close BCs: %i", bcRange.size());

    bool gA = true, gC = true;
    for (auto const& bc : bcRange) {
      if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits()))
        gA = false;
      if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits()))
        gC = false;
    }
    if (!gA && !gC)
      return 3;
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      return 4;
    }
    return gA && gC ? 2 : (gA ? 0 : 1);
  }
  template <typename TFwdTrack>
  int FwdTrkSelector(TFwdTrack const& fwdtrack)
  {
    if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3)
      return 1;
    else
      return 0;
  }

 private:
  TDatabasePDG* fPDG;
};

#endif // PWGUD_CORE_SGSELECTOR_H_
