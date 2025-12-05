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
//
/// \file   SGSelector.h
/// \author Alexander Bylinkin
/// \author Adam Matyja
/// \since  2023-11-21
/// \brief  Support class holding tools to work with Single Gap selection in central barrel UPCs
///

#ifndef PWGUD_CORE_SGSELECTOR_H_
#define PWGUD_CORE_SGSELECTOR_H_

#include "PWGUD/Core/SGCutParHolder.h"
#include "PWGUD/Core/UDHelpers.h"

#include "Common/CCDB/RCTSelectionFlags.h"

#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"

#include <cmath>
#include <memory>

template <typename BC>
struct SelectionResult {
  int value;    // The original integer return value
  std::shared_ptr<BC> bc; // Pointer to the BC object
};

namespace o2::aod::sgselector
{
enum TrueGap : int {
  NoGap = -1,        // no gap due to change of threshold(s) in any of FV0, FT0A, ZNA, FT0C, ZNC
  SingleGapA = 0,    // initially single gap at A side event
  SingleGapC = 1,    // initially single gap at C side event
  DoubleGap = 2,     // initially double gap event
  NoUpc = 3,         // initially no UPC event with default thresholds (FT0A=150, FT0C=50)
  TrkOutOfRange = 4, // to many tracks (>100 default)
  BadDoubleGap = 5   // unknows status of double gap check with changed thresholds
};
} // namespace o2::aod::sgselector

class SGSelector
{
 public:
  SGSelector() : myRCTChecker{"CBT"}, myRCTCheckerHadron{"CBT_hadronPID"}, myRCTCheckerZDC{"CBT", true}, myRCTCheckerHadronZDC{"CBT_hadronPID", true} {}

  template <typename CC, typename BCs, typename TCs, typename FWs>
  int Print(SGCutParHolder const& /*diffCuts*/, CC const& collision, BCs const& /*bcRange*/, TCs const& /*tracks*/, FWs const& /*fwdtracks*/)
  {
    LOGF(info, "Size of array %i", collision.size());
    return 1;
  }

  template <typename CC, typename BCs, typename BC>
  SelectionResult<BC> IsSelected(SGCutParHolder const& diffCuts, CC const& collision, BCs const& bcRange, BC const& oldbc)
  {
    //        LOGF(info, "Collision %f", collision.collisionTime());
    //        LOGF(info, "Number of close BCs: %i", bcRange.size());
    SelectionResult<BC> result;
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      result.value = o2::aod::sgselector::TrkOutOfRange; // 4
      result.bc = std::make_shared<BC>(oldbc);
      return result;
    }
    auto newbc = oldbc;
    auto newdgabc = oldbc;
    auto newdgcbc = oldbc;
    float tempampa = 0;
    float tempampc = 0;
    float ampc = 0;
    float ampa = 0;
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
    } // end of loop over bc range
    if (!gA && !gC) {
      result.value = o2::aod::sgselector::NoUpc; // gap = 3
      result.bc = std::make_shared<BC>(oldbc);
      return result;
    }
    if (gA && gC) { // loop once again for so-called DG events to get the most active FT0 BC
      for (auto const& bc : bcRange) {
        if (bc.has_foundFT0()) {
          tempampa = udhelpers::FT0AmplitudeA(bc.foundFT0());
          tempampc = udhelpers::FT0AmplitudeC(bc.foundFT0());
          if (tempampa > ampa) {
            ampa = tempampa;
            newdgabc = bc;
          }
          if (tempampc > ampc) {
            ampc = tempampc;
            newdgcbc = bc;
          }
        }
      }
      if (newdgabc != newdgcbc) {
        if (ampc / diffCuts.FITAmpLimits()[2] > ampa / diffCuts.FITAmpLimits()[1])
          newdgabc = newdgcbc;
      }
      newbc = newdgabc;
    }
    // LOGF(info, "Old BC: %i, New BC: %i",oldbc.globalBC(), newbc.globalBC());
    result.bc = std::make_shared<BC>(newbc);
    // result.value = gA && gC ? 2 : (gA ? 0 : 1);
    result.value = gA && gC ? o2::aod::sgselector::DoubleGap : (gA ? o2::aod::sgselector::SingleGapA : o2::aod::sgselector::SingleGapC);
    return result;
  }
  template <typename TFwdTrack>
  int FwdTrkSelector(TFwdTrack const& fwdtrack)
  {
    // if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3)
    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack || fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
      return 1;
    else
      return 0;
  }

  template <typename CC>
  int trueGap(CC const& collision, const float fv0, const float ft0a, const float ft0c, const float zdc_cut)
  {
    const float fit_cut[3] = {fv0, ft0a, ft0c};
    int gap = collision.gapSide();
    int true_gap = gap;
    // float FV0A, FT0A, FT0C, ZNA, ZNC;
    const float FV0A = collision.totalFV0AmplitudeA();
    const float FT0A = collision.totalFT0AmplitudeA();
    const float FT0C = collision.totalFT0AmplitudeC();
    const float ZNA = collision.energyCommonZNA();
    const float ZNC = collision.energyCommonZNC();
    if (gap == o2::aod::sgselector::SingleGapA) { // gap == 0
      if (FV0A > fit_cut[0] || FT0A > fit_cut[1] || ZNA > zdc_cut)
        true_gap = o2::aod::sgselector::NoGap;           // -1
    } else if (gap == o2::aod::sgselector::SingleGapC) { // gap == 1
      if (FT0C > fit_cut[2] || ZNC > zdc_cut)
        true_gap = o2::aod::sgselector::NoGap;          // -1
    } else if (gap == o2::aod::sgselector::DoubleGap) { // gap == 2
      if ((FV0A > fit_cut[0] || FT0A > fit_cut[1] || ZNA > zdc_cut) && (FT0C > fit_cut[2] || ZNC > zdc_cut)) {
        true_gap = o2::aod::sgselector::NoGap; // -1
      } else if ((FV0A > fit_cut[0] || FT0A > fit_cut[1] || ZNA > zdc_cut) && (FT0C <= fit_cut[2] && ZNC <= zdc_cut)) {
        true_gap = o2::aod::sgselector::SingleGapC; // 1
      } else if ((FV0A <= fit_cut[0] && FT0A <= fit_cut[1] && ZNA <= zdc_cut) && (FT0C > fit_cut[2] || ZNC > zdc_cut)) {
        true_gap = o2::aod::sgselector::SingleGapA; // 0
      } else if (FV0A <= fit_cut[0] && FT0A <= fit_cut[1] && ZNA <= zdc_cut && FT0C <= fit_cut[2] && ZNC <= zdc_cut) {
        true_gap = o2::aod::sgselector::DoubleGap; // 2
      } else {
        LOGF(info, "Something wrong with DG");
        true_gap = o2::aod::sgselector::BadDoubleGap; // 5
      }
    }
    return true_gap;
  }

  // check CBT flags
  template <typename CC>
  bool isCBTOk(CC const& collision)
  {
    if (myRCTChecker(collision))
      return true;
    return false;
  }

  // check CBT+hadronPID flags
  template <typename CC>
  bool isCBTHadronOk(CC const& collision)
  {
    if (myRCTCheckerHadron(collision))
      return true;
    return false;
  }

  // check CBT+ZDC flags
  template <typename CC>
  bool isCBTZdcOk(CC const& collision)
  {
    if (myRCTCheckerZDC(collision))
      return true;
    return false;
  }

  // check CBT+hadronPID+ZDC flags
  template <typename CC>
  bool isCBTHadronZdcOk(CC const& collision)
  {
    if (myRCTCheckerHadronZDC(collision))
      return true;
    return false;
  }

 private:
  o2::aod::rctsel::RCTFlagsChecker myRCTChecker;
  o2::aod::rctsel::RCTFlagsChecker myRCTCheckerHadron;
  o2::aod::rctsel::RCTFlagsChecker myRCTCheckerZDC;
  o2::aod::rctsel::RCTFlagsChecker myRCTCheckerHadronZDC;
};

#endif // PWGUD_CORE_SGSELECTOR_H_
