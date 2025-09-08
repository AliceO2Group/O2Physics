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
/// \file   CollisionTypeHelper.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Utility to handle the collision type from the GRP information
///

#ifndef COMMON_CORE_COLLISIONTYPEHELPER_H_
#define COMMON_CORE_COLLISIONTYPEHELPER_H_

#include <DataFormatsParameters/GRPLHCIFData.h>

#include <string>

namespace o2::common::core
{
// Container for the collision system type
struct CollisionSystemType {
  // Enum type for the collision system
  typedef int collType;

  static constexpr collType kCollSysUndef = -1; // Undefined collision system
  static constexpr collType kCollSyspp = 0;     // pp
  static constexpr collType kCollSysPbPb = 1;   // PbPb
  static constexpr collType kCollSysXeXe = 2;   // XeXe
  static constexpr collType kCollSyspPb = 3;    // pPb
  static constexpr collType kCollSysOO = 4;     // OO (Oxygen-Oxygen)
  static constexpr collType kCollSyspO = 5;     // pO (proton-Oxygen)
  static constexpr collType kCollSysNeNe = 6;   // NeNe (Neon-Neon)
  static constexpr collType kNCollSys = 7;      // Number of collision systems

  static std::string getCollisionSystemName(collType collSys);

  static int getCollisionTypeFromGrp(o2::parameters::GRPLHCIFData* grplhcif);
};

} // namespace o2::common::core

using CollisionSystemType = o2::common::core::CollisionSystemType;

#endif // COMMON_CORE_COLLISIONTYPEHELPER_H_
