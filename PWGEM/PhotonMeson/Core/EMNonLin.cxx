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

/// \file EMNonLin.cxx
/// \brief Source file of NonLin class for photons.
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "EMNonLin.h"

#include <cmath>

using namespace o2::pwgem::nonlin;

float EMNonLin::getCorrectionFactor(float inputCalibValue, PhotonType photonType, float cent)
{
  float param0 = 0, param1 = 0, param2 = 0, val = 1.f;
  switch (photonType) {
    case PhotonType::kEMC:
      if (cent >= 30 && cent <= 40) {
        param0 = -5.33426e-01f;
        param1 = 1.40144e-02f;
        param2 = -5.24434e-01f;
      } else {
        param0 = 0.f;
        param1 = 0.f;
        param2 = 0.f;
      }
      break;
    case PhotonType::kPCM:
      if (cent >= 0 && cent <= 100) {
        param0 = 10.7203f;
        param1 = 0.0383968f;
        param2 = 10.6025f;
      } else {
        param0 = 0.f;
        param1 = 0.f;
        param2 = 0.f;
      }
      break;
    case PhotonType::kPHOS:
      param0 = 0.f;
      param1 = 0.f;
      param2 = 0.f;
      break;
  }

  val = (1.f + param0 / inputCalibValue + param1 / std::pow(inputCalibValue, 2.f)) / (1.f + param2 / inputCalibValue);
  return val;
}
