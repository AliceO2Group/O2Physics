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
DetLayer::DetLayer(const TString& name,
                   float r,
                   float z,
                   float x0,
                   float xrho,
                   float resRPhi,
                   float resZ,
                   float eff,
                   int type)
  : mName(name),
    mR(r),
    mZ(z),
    mX0(x0),
    mXrho(xrho),
    mResRPhi(resRPhi),
    mResZ(resZ),
    mEff(eff),
    mType(type)
{
}

DetLayer::DetLayer(const DetLayer& other)
  : mName(other.mName), mR(other.mR), mZ(other.mZ), mX0(other.mX0), mXrho(other.mXrho), mResRPhi(other.mResRPhi), mResZ(other.mResZ), mEff(other.mEff), mType(other.mType)
{
}

void DetLayer::addDeadPhiRegion(float phiStart, float phiEnd)
{
  static constexpr float DefaultValue = 2.f;
  static constexpr float PhiTolerance = 1e-4f;
  if (mDeadPhiRegions == nullptr) {
    mDeadPhiRegions = new TGraph();
    mDeadPhiRegions->SetNameTitle(Form("deadPhiRegions_%s", mName.Data()), Form("Dead phi regions for layer %s", mName.Data()));
    mDeadPhiRegions->AddPoint(0, DefaultValue);
    mDeadPhiRegions->AddPoint(o2::constants::math::TwoPI, DefaultValue);
  }
  if (phiStart < 0 || phiStart >= o2::constants::math::TwoPI || phiEnd < 0 || phiEnd >= o2::constants::math::TwoPI) {
    LOG(fatal) << "Cannot add dead phi region with invalid range [" << phiStart << ", " << phiEnd << "] to layer " << mName;
    return;
  }
  mDeadPhiRegions->AddPoint(phiStart, DefaultValue);
  mDeadPhiRegions->AddPoint(phiEnd, DefaultValue);
  mDeadPhiRegions->AddPoint(phiStart + PhiTolerance, 0.f);
  mDeadPhiRegions->AddPoint(phiEnd - PhiTolerance, 0.f);
  mDeadPhiRegions->Sort();
}

void DetLayer::setDeadPhiRegions(TGraph* graph)
{
  LOG(debug) << "Setting dead phi regions for layer " << mName << " with graph " << (graph ? graph->GetName() : "nullptr");
  if (mDeadPhiRegions != nullptr) {
    LOG(warning) << "Overriding existing dead phi regions for layer " << mName;
    delete mDeadPhiRegions;
  }
  mDeadPhiRegions = graph;
  if (mDeadPhiRegions->GetN() == 0) {
    LOG(warning) << "Dead phi regions graph for layer " << mName << " is empty, clearing dead regions";
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
        LOG(fatal) << "Dead phi regions graph for layer " << mName << " has invalid x value " << x << " at point " << i << ", first point should be 0 and last point should be 2PI";
      }
      LOG(debug) << "Point " << i << ": (" << x << ", " << y << ")";
      if (x < 0 || x > o2::constants::math::TwoPI) {
        LOG(fatal) << "Dead phi regions graph for layer " << mName << " has invalid x value " << x << " at point " << i;
      }
      static constexpr float EffValidLowValue = 0.f;
      static constexpr float EffValidHighValue = 2.f;
      if (y != EffValidLowValue && y != EffValidHighValue) {
        LOG(fatal) << "Dead phi regions graph for layer " << mName << " has invalid y value " << y << " at point " << i << ", should be 0 or 2";
      }
    }
  } else {
    LOG(info) << "Cleared dead phi regions for layer " << mName;
  }
}

std::string DetLayer::toString() const
{
  std::string out = "";
  out.append("DetLayer: ");
  out.append(mName.Data());
  out.append(" | r: ");
  out.append(std::to_string(mR));
  out.append(" cm | z: ");
  out.append(std::to_string(mZ));
  out.append(" cm | x0: ");
  out.append(std::to_string(mX0));
  out.append(" cm | xrho: ");
  out.append(std::to_string(mXrho));
  out.append(" g/cm^3 | resRPhi: ");
  out.append(std::to_string(mResRPhi));
  out.append(" cm | resZ: ");
  out.append(std::to_string(mResZ));
  out.append(" cm | eff: ");
  out.append(std::to_string(mEff));
  out.append(" | type: ");
  switch (mType) {
    case kLayerInert:
      out.append("Inert");
      break;
    case kLayerSilicon:
      out.append("Silicon");
      break;
    case kLayerGas:
      out.append("Gas/TPC");
      break;
    case kLayerTOF:
      out.append("TOF");
      break;
    case kLayerVertex:
      out.append("Vertex");
      break;
    default:
      out.append("Unknown");
      break;
  }
  return out;
}

} // namespace o2::fastsim
