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
/// \file   TOFReso.h
/// \author Nicolo' Jacazio
/// \since  07/08/2020
/// \brief  Implementation for the TOF PID response of the expected times resolution
///

#ifndef O2_ANALYSIS_PID_TOFRESO_H_
#define O2_ANALYSIS_PID_TOFRESO_H_

// O2 includes
#include "ReconstructionDataFormats/PID.h"

namespace o2::pid::tof
{

class TOFReso : public Parametrization
{
 public:
  TOFReso() : Parametrization("TOFReso", 5){};
  ~TOFReso() override = default;
  /// Operator to compute the expected value of the TOF Resolution
  /// \param x Array with the input used to compute the response:
  /// x[0] -> track momentum
  /// x[1] -> TOF signal
  /// x[2] -> event time resolution
  /// x[3] -> particle mass
  float operator()(const float* x) const override
  {
    const float mom = abs(x[0]);
    if (mom <= 0) {
      return -999;
    }
    const float time = x[1];
    const float evtimereso = x[2];
    const float mass = x[3];
    const float dpp = mParameters[0] + mParameters[1] * mom + mParameters[2] * mass / mom; // mean relative pt resolution;
    const float sigma = dpp * time / (1. + mom * mom / (mass * mass));
    return sqrt(sigma * sigma + mParameters[3] * mParameters[3] / mom / mom + mParameters[4] * mParameters[4] + evtimereso * evtimereso);
  }
  ClassDefOverride(TOFReso, 1);
};

} // namespace o2::pid::tof

#endif
