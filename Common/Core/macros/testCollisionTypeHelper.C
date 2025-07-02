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

/// \file   testCollisionTypeHelper.C
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Test the CollisionTypeHelper functionality

#include "Common/Core/CollisionTypeHelper.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

void testCollisionTypeHelper(int runNumber = 544124)
{

  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  ccdb.setURL("http://alice-ccdb.cern.ch");
  o2::parameters::GRPLHCIFData* grpo = ccdb.getSpecificForRun<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF",
                                                                                            runNumber);
  grpo->print();
  int collsys = o2::common::core::CollisionSystemType::getCollisionTypeFromGrp(grpo);
}
