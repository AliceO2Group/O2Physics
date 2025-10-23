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
/// \file   DetLayer.cxx
/// \author David Dobrigkeit Chinellato
/// \since  11/03/2021
/// \brief  Basic struct to hold information regarding a detector layer to be used in fast simulation
///

#include "DetLayer.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>

#include <string>
#include <vector>

namespace o2::fastsim
{

// Parametric constructor
DetLayer::DetLayer(const TString& name_,
                   float r_,
                   float z_,
                   float x0_,
                   float xrho_,
                   float resRPhi_,
                   float resZ_,
                   float eff_,
                   int type_)
  : name(name_),
    r(r_),
    z(z_),
    x0(x0_),
    xrho(xrho_),
    resRPhi(resRPhi_),
    resZ(resZ_),
    eff(eff_),
    type(type_)
{
}

DetLayer::DetLayer(const DetLayer& other)
  : name(other.name), r(other.r), z(other.z), x0(other.x0), xrho(other.xrho), resRPhi(other.resRPhi), resZ(other.resZ), eff(other.eff), type(other.type)
{
}

void DetLayer::addDeadPhiRegion(float phiStart, float phiEnd)
{
  static constexpr float kDefaultValue = 2.f;
  static constexpr float kPhiTolerance = 1e-4f;
  if (mDeadPhiRegions == nullptr) {
    mDeadPhiRegions = new TGraph();
    mDeadPhiRegions->SetNameTitle(Form("deadPhiRegions_%s", name.Data()), Form("Dead phi regions for layer %s", name.Data()));
    mDeadPhiRegions->AddPoint(0, kDefaultValue);
    mDeadPhiRegions->AddPoint(o2::constants::math::TwoPI, kDefaultValue);
  }
  if (phiStart < 0 || phiStart >= o2::constants::math::TwoPI || phiEnd < 0 || phiEnd >= o2::constants::math::TwoPI) {
    LOG(fatal) << "Cannot add dead phi region with invalid range [" << phiStart << ", " << phiEnd << "] to layer " << name;
    return;
  }
  mDeadPhiRegions->AddPoint(phiStart, kDefaultValue);
  mDeadPhiRegions->AddPoint(phiEnd, kDefaultValue);
  mDeadPhiRegions->AddPoint(phiStart + kPhiTolerance, 0.f);
  mDeadPhiRegions->AddPoint(phiEnd - kPhiTolerance, 0.f);
  mDeadPhiRegions->Sort();
}

void DetLayer::setDeadPhiRegions(TGraph* graph)
{
  LOG(debug) << "Setting dead phi regions for layer " << name << " with graph " << (graph ? graph->GetName() : "nullptr");
  if (mDeadPhiRegions != nullptr) {
    LOG(warning) << "Overriding existing dead phi regions for layer " << name;
    delete mDeadPhiRegions;
  }
  mDeadPhiRegions = graph;
  if (mDeadPhiRegions->GetN() == 0) {
    LOG(warning) << "Dead phi regions graph for layer " << name << " is empty, clearing dead regions";
    mDeadPhiRegions = nullptr;
    return; // cleared the dead regions
  }
  // Check sanity of the graph
  if (mDeadPhiRegions != nullptr) {
    for (int i = 0; i < mDeadPhiRegions->GetN(); i++) {
      const float x = mDeadPhiRegions->GetX()[i];
      const float y = mDeadPhiRegions->GetY()[i];
      // First point has to be at 0, last point has to be at 2PI
      if ((i == 0 && x != 0.f) || (i == mDeadPhiRegions->GetN() - 1 && x != o2::constants::math::TwoPI)) {
        LOG(fatal) << "Dead phi regions graph for layer " << name << " has invalid x value " << x << " at point " << i << ", first point should be 0 and last point should be 2PI";
      }
      LOG(debug) << "Point " << i << ": (" << x << ", " << y << ")";
      if (x < 0 || x > o2::constants::math::TwoPI) {
        LOG(fatal) << "Dead phi regions graph for layer " << name << " has invalid x value " << x << " at point " << i;
      }
      if (y != 0.f && y != 2.f) {
        LOG(fatal) << "Dead phi regions graph for layer " << name << " has invalid y value " << y << " at point " << i << ", should be 0 or 2";
      }
    }
  } else {
    LOG(info) << "Cleared dead phi regions for layer " << name;
  }
}

std::string DetLayer::toString() const
{
  std::string out = "";
  out.append("DetLayer: ");
  out.append(name.Data());
  out.append(" | r: ");
  out.append(std::to_string(r));
  out.append(" cm | z: ");
  out.append(std::to_string(z));
  out.append(" cm | x0: ");
  out.append(std::to_string(x0));
  out.append(" cm | xrho: ");
  out.append(std::to_string(xrho));
  out.append(" g/cm^3 | resRPhi: ");
  out.append(std::to_string(resRPhi));
  out.append(" cm | resZ: ");
  out.append(std::to_string(resZ));
  out.append(" cm | eff: ");
  out.append(std::to_string(eff));
  out.append(" | type: ");
  switch (type) {
    case layerInert:
      out.append("Inert");
      break;
    case layerSilicon:
      out.append("Silicon");
      break;
    case layerGas:
      out.append("Gas/TPC");
      break;
    default:
      out.append("Unknown");
      break;
  }
  return out;
}

} // namespace o2::fastsim
