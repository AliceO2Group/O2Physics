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

#include <cmath>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGCutParHolder.h"

template <typename BC>
struct SelectionResult {
  int value; // The original integer return value
  BC* bc;    // Pointer to the BC object
};

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

  template <typename CC, typename BCs, typename BC>
  SelectionResult<BC> IsSelected(SGCutParHolder diffCuts, CC& collision, BCs& bcRange, BC& oldbc)
  {
    //        LOGF(info, "Collision %f", collision.collisionTime());
    //        LOGF(info, "Number of close BCs: %i", bcRange.size());
    SelectionResult<BC> result;
    result.bc = &oldbc;
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      result.value = 4;
      return result;
    }
    auto newbc = oldbc;
    bool gA = true, gC = true;
    for (auto const& bc : bcRange) {
      if (!udhelpers::cleanFITA(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        if (gA)
          newbc = bc;
        if (!gA && std::abs(static_cast<int64_t>(bc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newbc.globalBC() - oldbc.globalBC())))
          newbc = bc;
        gA = false;
      }
      if (!udhelpers::cleanFITC(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        if (gC)
          newbc = bc;
        if (!gC && std::abs(static_cast<int64_t>(bc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newbc.globalBC() - oldbc.globalBC())))
          newbc = bc;
        gC = false;
      }
    }
    result.bc = &newbc;
    if (!gA && !gC) {
      result.value = 3;
      return result;
    }
    // LOGF(info, "Old BC: %i, New BC: %i",oldbc.globalBC(), newbc.globalBC());
    result.value = gA && gC ? 2 : (gA ? 0 : 1);
    return result;
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
