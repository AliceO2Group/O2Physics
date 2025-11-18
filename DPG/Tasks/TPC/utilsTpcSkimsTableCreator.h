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

/// \file utilsTpcSkimsTableCreator.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Helper functions used both in V0 and TOF structs of tpcSkimsTableCreator.cxx

#ifndef DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
#define DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_

#include <TRandom3.h>

namespace o2::dpg_tpcskimstablecreator
{
enum {
  TrackSelectionNoCut = 0,
  TrackSelectionGlobalTrack,
  TrackSelectionTrackWoPtEta,
  TrackSelectionGlobalTrackWoDCA,
  TrackSelectionQualityTracks,
  TrackSelectionInAcceptanceTracks
};

enum {
  EventSelectionNo = 0,
  EventSelectionRun2,
  EventSelectionRun3
};

enum {
  ModeStandard = 0,
  ModeWithdEdxTrkQA,
  ModeWithTrkQA
};

/// Event selection
template <typename CollisionType>
inline bool isEventSelected(const CollisionType& collision, const int applyEvSel)
{
  if ((applyEvSel == EventSelectionRun2 && !collision.sel7()) || (applyEvSel == EventSelectionRun3 && !collision.sel8())) {
    return false;
  } else {
    return true;
  }
};

/// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV)
/// as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
inline bool downsampleTsalisCharged(TRandom3* fRndm, const double pt, const double factor1Pt, const double mass, const double sqrtSNN, const double maxPt = 1e9)
{
  if (factor1Pt < 0. || pt > maxPt) {
    return true;
  }

  auto tsalisCharged = [&](const double pT) {
    const double a = 6.81, b = 59.24;
    const double c = 0.082, d = 0.151;
    const double mt = std::sqrt(mass * mass + pT * pT);
    const double n = a + b / sqrtSNN;
    const double t = c + d / sqrtSNN;
    const double p0 = n * t;
    const double result = std::pow((1. + mt / p0), -n);
    return result;
  };

  const double prob = tsalisCharged(pt) * pt;
  const double probNorm = tsalisCharged(1.);
  if ((fRndm->Rndm() * ((prob / probNorm) * pt * pt)) > factor1Pt) {
    return false;
  } else {
    return true;
  }
};

// Track selection
template <typename TrackType>
inline bool isTrackSelected(const TrackType& track, const int trackSelection)
{
  bool isSelected{false};
  isSelected |= trackSelection == TrackSelectionNoCut;
  isSelected |= (trackSelection == TrackSelectionGlobalTrack) && track.isGlobalTrack();
  isSelected |= (trackSelection == TrackSelectionTrackWoPtEta) && track.isGlobalTrackWoPtEta();
  isSelected |= (trackSelection == TrackSelectionGlobalTrackWoDCA) && track.isGlobalTrackWoDCA();
  isSelected |= (trackSelection == TrackSelectionQualityTracks) && track.isQualityTrack();
  isSelected |= (trackSelection == TrackSelectionInAcceptanceTracks) && track.isInAcceptanceTrack();

  return isSelected;
}

} // namespace o2::dpg_tpcskimstablecreator
#endif // DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
