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

/// \commonly used for PCM analyses.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_

#include "Framework/AnalysisTask.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

//_______________________________________________________________________
TrackSelection CreateTPCOnlyTrackCuts(float max_eta, int min_tpc_ncr, float max_chi2_tpc)
{
  TrackSelection selectedTracks;
  selectedTracks.SetPtRange(0.01f, 1e10f);
  selectedTracks.SetEtaRange(-max_eta, max_eta);
  selectedTracks.SetMinNCrossedRowsTPC(min_tpc_ncr);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.6f);
  selectedTracks.SetMaxChi2PerClusterTPC(max_chi2_tpc);
  selectedTracks.SetMaxDcaXY(2.4f);
  // selectedTracks.SetMaxDcaZ(3.2f);
  return selectedTracks;
}
//_______________________________________________________________________
bool checkAP(float alpha, float qt)
{
  const float alpha_max = 0.95;
  const float qt_max = 0.05;
  float ellipse = pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2);
  if (ellipse < 1.0) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
//_______________________________________________________________________
//_______________________________________________________________________
#endif // PWGEM_PHOTONMESON_UTILS_PCMUTILITIES_H_
