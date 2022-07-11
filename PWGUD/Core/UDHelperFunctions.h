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
///
/// \brief
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#ifndef O2_ANALYSISUDHEPLER_H_
#define O2_ANALYSISUDHEPLER_H_

#include "Framework/Logger.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;

// .............................................................................
// return net charge of PV tracks
template <typename TCs>
int8_t netCharge(TCs tracks)
{
  int8_t nch = 0;
  for (auto track : tracks) {
    if (track.isPVContributor()) {
      nch += track.sign();
    }
  }
  return nch;
}

// .............................................................................
// return fraction of PV tracks with a TOF hit
template <typename TCs>
float rPVtrwTOF(TCs tracks, int nPVTracks)
{
  float rpvrwTOF = 0.;
  for (auto& track : tracks) {
    if (track.isPVContributor() && track.hasTOF()) {
      rpvrwTOF += 1.;
    }
  }
  if (nPVTracks > 0) {
    rpvrwTOF /= nPVTracks;
  }
  return rpvrwTOF;
}

// -----------------------------------------------------------------------------
#endif // O2_ANALYSISUDHEPLER_H_
