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

/// \file   testFastTracker.C
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Test the FastTracker functionality

#include "ALICE3/Core/FastTracker.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPLHCIFData.h>

#include <fairlogger/Logger.h>

void testFastTracker(std::string geometryFile = "a3geo.ini")
{

  fair::Logger::SetConsoleSeverity(fair::Severity::debug);

  // auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  // ccdb.setURL("http://alice-ccdb.cern.ch");
  o2::fastsim::FastTracker fastTracker;
  fastTracker.AddGenericDetector(geometryFile);
  // fastTracker.AddGenericDetector(geometryFile, &ccdb);
  fastTracker.Print();
}
