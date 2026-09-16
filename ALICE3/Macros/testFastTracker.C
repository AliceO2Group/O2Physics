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
#include "ALICE3/Core/GeometryContainer.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/Logger.h>

#include <string>

void testFastTracker(const std::string& geometryFile = "Configuration/a3geo.ini")
{

  fair::Logger::SetConsoleSeverity(fair::Severity::debug);

  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  ccdb.setURL("http://alice-ccdb.cern.ch");

  o2::fastsim::GeometryContainer geometryContainer;
  geometryContainer.setCcdbManager(&ccdb);
  geometryContainer.addEntry(geometryFile);
  const int nGeometries = geometryContainer.getNumberOfConfigurations();
  const float magneticField = geometryContainer.getFloatValue(0, "global", "magneticfield");

  for (int icfg = 0; icfg < nGeometries; ++icfg) {
    const o2::fastsim::GeometryEntry& geometry = geometryContainer.getEntry(icfg);
    o2::fastsim::FastTracker fastTracker;

    fastTracker.setMagneticField(magneticField);
    fastTracker.addGenericDetector(geometry, &ccdb);

    fastTracker.print();
  }
}
