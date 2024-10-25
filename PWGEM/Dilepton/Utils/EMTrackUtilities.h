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
template <typename T>
float sigmaDCA3D(T const& track)
{
  return 1. / dca3DinSigma(track) * std::sqrt(pow(track.dcaXY(), 2) + pow(track.dcaZ(), 2));

  // float cYY = track.cYY();
  // float cZZ = track.cZZ();
  // float cZY = track.cZY();
  // float dcaXY = track.dcaXY();                                                                                                                                // in cm
  // float dcaZ = track.dcaZ();                                                                                                                                  // in cm
  // return 1.f / std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ) * std::sqrt(pow(dcaXY * std::sqrt(cYY), 2) + pow(dcaZ * std::sqrt(cZZ), 2) + 2.f * dcaXY * dcaZ * cZY); // DCA 3D resolution at midrapidity
}
//_______________________________________________________________________
template <typename T>
float sigmaDCAXY(T const& track)
{
  return 1. / fwdDcaXYinSigma(track) * std::sqrt(pow(track.fwdDcaX(), 2) + pow(track.fwdDcaY(), 2));

  // float cXX = track.cXX();
  // float cYY = track.cYY();
  // float cXY = track.cXY();
  // float dcaX = track.fwdDcaX();                                                                                                                           // in cm
  // float dcaY = track.fwdDcaY();                                                                                                                           // in cm
  // return 1.f / std::sqrt(dcaX * dcaX + dcaY * dcaY) * std::sqrt(pow(dcaX * std::sqrt(cXX), 2) + pow(dcaY * std::sqrt(cYY), 2) + 2.f * dcaX * dcaY * cXY); // DCA XY resolution at forward rapidity
}
//_______________________________________________________________________
template <typename T>
float sigmaPt(T const& track)
{
  return std::sqrt(track.c1Pt21Pt2()) / std::pow(track.signed1Pt(), 2); // pT resolution
}
//_______________________________________________________________________
template <typename T>
float sigmaPhi(T const& track)
{
  return std::sqrt(track.cSnpSnp()) / std::sqrt(1.f - std::pow(track.snp(), 2)); // phi resolution
}
//_______________________________________________________________________
template <typename T>
float sigmaTheta(T const& track)
{
  return std::sqrt(track.cTglTgl()) / (1.f + std::pow(track.tgl(), 2)); // theta resolution = lambda resolution. // lambda = pi/2 - theta. theta is polar angle.
}
//_______________________________________________________________________
template <typename T>
float sigmaEta(T const& track)
{
  return std::sqrt(track.cTglTgl()) / std::sqrt(1.f + std::pow(track.tgl(), 2));
}
//_______________________________________________________________________
template <typename T>
float sigmaP(T const& track)
{
  // p = 1/1/pT x 1/cos(lambda);
  return std::sqrt(std::pow(1.f / track.signed1Pt(), 4) * ((1.f + std::pow(track.tgl(), 2)) * track.c1Pt21Pt2() + 1.f / (1.f + std::pow(track.tgl(), 2)) * std::pow(track.signed1Pt() * track.tgl(), 2) * track.cTglTgl() - 2.f * track.signed1Pt() * track.tgl() * track.c1PtTgl()));
}
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::emtrackutil
#endif // PWGEM_DILEPTON_UTILS_EMTRACKUTILITIES_H_
