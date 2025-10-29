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
/// \file   CollisionTypeHelper.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Utility to handle the collision type from the GRP information
///

#include "Common/Core/CollisionTypeHelper.h"

#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/Logger.h>

#include <string>

std::string o2::common::core::CollisionSystemType::getCollisionSystemName(collType collSys)
{
  switch (collSys) {
    case kCollSyspp:
      return "pp";
    case kCollSysPbPb:
      return "PbPb";
    case kCollSysXeXe:
      return "XeXe";
    case kCollSyspPb:
      return "pPb";
    case kCollSysOO:
      return "OO";
    case kCollSyspO:
      return "pO";
    case kCollSysNeNe:
      return "NeNe";
    case kCollSysUndef:
      return "Undefined";
    default:
      LOG(warning) << "Undefined collision system type: " << collSys;
      return "Undefined";
  }
}

int o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(o2::parameters::GRPLHCIFData* grplhcif)
{
  const int zBeamA = grplhcif->getBeamZ(o2::constants::lhc::BeamDirection::BeamA);
  const int zBeamC = grplhcif->getBeamZ(o2::constants::lhc::BeamDirection::BeamC);
  LOG(debug) << "Collision system Z: " << zBeamA << " * " << zBeamC << " detected = " << zBeamA * zBeamC;
  switch (zBeamA * zBeamC) {
    case 1: // pp 1*1
      return kCollSyspp;
    case 6724: // Pb-Pb 82*82
      return kCollSysPbPb;
    case 225: // Xe-Xe 54*54
      return kCollSysXeXe;
    case 82: // p-Pb 82*1
      return kCollSyspPb;
    case 64: // O-O 8*8
      return kCollSysOO;
    case 8: // p-O 8*1
      return kCollSyspO;
    case 100: // Ne-Ne 10*10
      return kCollSysNeNe;
    default:
      LOG(fatal) << "Undefined collision system in getCollisionTypeFromGrp with Z of BeamA = " << zBeamA << " and Z of BeamC = " << zBeamC;
      return kCollSysUndef;
  }
  return kCollSysUndef;
}
