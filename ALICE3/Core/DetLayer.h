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
/// \file   DetLayer.h
/// \author David Dobrigkeit Chinellato
/// \since  11/03/2021
/// \brief  Basic struct to hold information regarding a detector layer to be used in fast simulation
///

#ifndef ALICE3_CORE_DETLAYER_H_
#define ALICE3_CORE_DETLAYER_H_

#include "TString.h"

namespace o2::fastsim
{

struct DetLayer {
  // TString for holding name
  TString name;

  // position variables
  float r; // radius in centimeters
  float z; // z dimension in centimeters

  // material variables
  float x0;   // radiation length
  float xrho; // density

  // resolution variables for active layers
  float resRPhi; // RPhi resolution in centimeters
  float resZ;    // Z resolution in centimeters

  // efficiency
  float eff; // detection efficiency

  // layer type
  int type; // 0: undefined/inert, 1: silicon, 2: gas/tpc
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_DETLAYER_H_
