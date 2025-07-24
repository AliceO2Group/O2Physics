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
/// \file   TOFResoALICE3.h
/// \author Nicolo' Jacazio
/// \since  11/03/2021
/// \brief  Implementation for the TOF PID response of the expected times resolution
///

#include "ALICE3/Core/TOFResoALICE3.h"

namespace o2::pid::tof
{

float TOFResoALICE3Param(const float& momentum, const float& momentumError, const float& evtimereso, const float& length, const float& mass, const Parameters& parameters)
{
  if (momentum <= 0) {
    return -999.f;
  }

  const float p2 = momentum * momentum;
  const float Lc = length / 0.0299792458f;
  const float mass2 = mass * mass;
  const float ep = momentumError * momentum;
  // const float ep = momentumError * p2;
  const float etexp = Lc * mass2 / p2 / sqrt(mass2 + p2) * ep;
  return sqrt(etexp * etexp + parameters[0] * parameters[0] + evtimereso * evtimereso);
}

} // namespace o2::pid::tof
