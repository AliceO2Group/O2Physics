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

#include <vector>
#include <string>

#include "DetLayer.h"

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
