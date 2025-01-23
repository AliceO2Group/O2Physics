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
/// \file   PIDTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  02/07/2020
/// \brief  Implementation of the TOF detector response for PID
///

#include "PIDTOF.h"
#include <string>

namespace o2::pid::tof
{

void TOFResoParamsV3::setResolutionParametrizationRun2(std::unordered_map<std::string, float> const& pars)
{
  std::array<std::string, 13> paramNames{"TrkRes.Pi.P0", "TrkRes.Pi.P1", "TrkRes.Pi.P2", "TrkRes.Pi.P3", "time_resolution",
                                         "TrkRes.Ka.P0", "TrkRes.Ka.P1", "TrkRes.Ka.P2", "TrkRes.Ka.P3",
                                         "TrkRes.Pr.P0", "TrkRes.Pr.P1", "TrkRes.Pr.P2", "TrkRes.Pr.P3"};
  // Now we override the parametrization to use the Run 2 one
  for (int i = 0; i < 9; i++) {
    if (mResolution[i]) {
      delete mResolution[i];
    }
    mResolution[i] = new TF2(Form("tofResTrack.%s_Run2", particleNames[i]), "-10", 0., 20, -1, 1.); // With negative values the old one is used
  }
  // Print the map
  for (const auto& [key, value] : pars) {
    LOG(info) << "Key: " << key << " Value: " << value;
  }
  for (int i = 0; i < 13; i++) {
    setParameter(i, pars.at(paramNames[i]));
  }
}

} // namespace o2::pid::tof
