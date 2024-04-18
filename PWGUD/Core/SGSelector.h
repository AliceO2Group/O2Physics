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
    auto newznabc = oldbc;
    auto newzncbc = oldbc;
    bool gA = true, gC = true;
    bool gzA = true, gzC = true;
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
    if (!gA && !gC) {
      result.value = 3;
      return result;
    }
    for (auto const& bc : bcRange) {
      if (bc.has_zdc()) {
        auto zdc = bc.zdc();
        if (std::abs(static_cast<float>(zdc.timeZNA())) < 2 && zdc.energyCommonZNA() > 0) {
          if (gzA) {
            newznabc = bc;
          }
          if (!gzA && std::abs(static_cast<int64_t>(bc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newznabc.globalBC() - oldbc.globalBC()))) {
            newznabc = bc;
          }
          gzA = false;
        }
        if (std::abs(static_cast<float>(zdc.timeZNC())) < 2 && zdc.energyCommonZNC() > 0) {
          if (gzC) {
            newzncbc = bc;
          }
          if (!gzC && std::abs(static_cast<int64_t>(bc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newzncbc.globalBC() - oldbc.globalBC()))) {
            newzncbc = bc;
          }
          gzC = false;
        }
      }
    }
    if (gA && gC) {
      if (!gzA && !gzC) {
        if (std::abs(static_cast<int64_t>(newznabc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newzncbc.globalBC() - oldbc.globalBC()))) {
          newzncbc = newznabc;
          newbc = newznabc;
        } else {
          newznabc = newzncbc;
          newbc = newzncbc;
        }
      } else if (!gzA) {
        newzncbc = newznabc;
        newbc = newznabc;
      } else if (!gzC) {
        newznabc = newzncbc;
        newbc = newzncbc;
      }
    } else if (!gA) {
      if (!gzA) {
        if (newbc.globalBC() == newznabc.globalBC()) {
          newzncbc = newznabc;
        } else if (std::abs(static_cast<int64_t>(newznabc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newbc.globalBC() - oldbc.globalBC()))) {
          newzncbc = newznabc;
          newbc = newznabc;
        } else {
          newzncbc = newbc;
          newznabc = newbc;
        }
      } else {
        newzncbc = newbc;
        newznabc = newbc;
      }
    } else if (!gC) {
      if (!gzC) {
        if (newbc.globalBC() == newzncbc.globalBC()) {
          newznabc = newzncbc;
        } else if (std::abs(static_cast<int64_t>(newzncbc.globalBC() - oldbc.globalBC())) < std::abs(static_cast<int64_t>(newbc.globalBC() - oldbc.globalBC()))) {
          newznabc = newzncbc;
          newbc = newzncbc;
        } else {
          newzncbc = newbc;
          newznabc = newbc;
        }

      } else {
        newzncbc = newbc;
        newznabc = newbc;
      }
    }
    // LOGF(info, "Old BC: %i, New BC: %i",oldbc.globalBC(), newbc.globalBC());
    result.bc = &newbc;
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

  template <typename CC>
  int trueGap(CC& collision, float fv0_cut, float zdc_cut)
  {
    int gap = collision.gapSide();
    int true_gap = gap;
    if (gap == 0) {
      if (collision.totalFV0AmplitudeA() > fv0_cut || collision.energyCommonZNA() > zdc_cut)
        true_gap = -1;
    } else if (gap == 1) {
      if (collision.energyCommonZNC() > zdc_cut)
        true_gap = -1;
    } else if (gap == 2) {
      if ((collision.totalFV0AmplitudeA() > fv0_cut || collision.energyCommonZNA() > zdc_cut) && (collision.energyCommonZNC() > zdc_cut))
        true_gap = -1;
      else if (collision.totalFV0AmplitudeA() > fv0_cut || collision.energyCommonZNA() > zdc_cut)
        true_gap = 1;
      else if (collision.energyCommonZNC() > zdc_cut)
        true_gap = 0;
    }
    return true_gap;
  }

 private:
  TDatabasePDG* fPDG;
};

#endif // PWGUD_CORE_SGSELECTOR_H_
