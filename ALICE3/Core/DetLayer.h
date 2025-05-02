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

#include <string>

#include "TString.h"

namespace o2::fastsim
{

struct DetLayer {
 public:
  // Default constructor
  DetLayer() = default;
  // Parametric constructor
  DetLayer(const TString& name_, float r_, float z_, float x0_, float xrho_,
           float resRPhi_ = 0.0f, float resZ_ = 0.0f, float eff_ = 0.0f, int type_ = layerInert);
  // Copy constructor
  DetLayer(const DetLayer& other);

  // Setters
  void setName(const TString& name_) { name = name_; }
  void setRadius(float r_) { r = r_; }
  void setZ(float z_) { z = z_; }
  void setRadiationLength(float x0_) { x0 = x0_; }
  void setDensity(float xrho_) { xrho = xrho_; }
  void setResolutionRPhi(float resRPhi_) { resRPhi = resRPhi_; }
  void setResolutionZ(float resZ_) { resZ = resZ_; }
  void setEfficiency(float eff_) { eff = eff_; }
  void setType(int type_) { type = type_; }

  // Getters
  float getRadius() const { return r; }
  float getZ() const { return z; }
  float getRadiationLength() const { return x0; }
  float getDensity() const { return xrho; }
  float getResolutionRPhi() const { return resRPhi; }
  float getResolutionZ() const { return resZ; }
  float getEfficiency() const { return eff; }
  int getType() const { return type; }
  const TString& getName() const { return name; }

  // Check layer type
  bool isInert() const { return type == layerInert; }
  bool isSilicon() const { return type == layerSilicon; }
  bool isGas() const { return type == layerGas; }

  // Utilities
  std::string toString() const;
  friend std::ostream& operator<<(std::ostream& os, const DetLayer& layer)
  {
    os << layer.toString();
    return os;
  }

 private:
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
  int type;                              // 0: undefined/inert, 1: silicon, 2: gas/tpc
  static constexpr int layerInert = 0;   // inert/undefined layer
  static constexpr int layerSilicon = 1; // silicon layer
  static constexpr int layerGas = 2;     // gas/tpc layer
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_DETLAYER_H_
