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

/// \commonly used to calculate track variables
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_
#define PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_

#include <string>
#include <vector>
#include <algorithm>

//_______________________________________________________________________
namespace o2::aod::pwgem::dilepton::utils::emtrackutil
{
template <typename T>
float dca3DinSigma(T const& track)
{
  float cYY = track.cYY();
  float cZZ = track.cZZ();
  float cZY = track.cZY();
  float dcaXY = track.dcaXY(); // in cm
  float dcaZ = track.dcaZ();   // in cm

  float det = cYY * cZZ - cZY * cZY; // determinant
  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::abs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
  }
}
//_______________________________________________________________________
template <typename T>
float fwdDcaXYinSigma(T const& track)
{
  float cXX = track.cXX();
  float cYY = track.cYY();
  float cXY = track.cXY();
  float dcaX = track.fwdDcaX(); // in cm
  float dcaY = track.fwdDcaY(); // in cm

  float det = cXX * cYY - cXY * cXY; // determinant
  if (det < 0) {
    return 999.f;
  } else {
    return std::sqrt(std::abs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2. * dcaX * dcaY * cXY) / det / 2.)); // dca xy in sigma
  }
}
//_______________________________________________________________________
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::emtrackutil
#endif // PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_
