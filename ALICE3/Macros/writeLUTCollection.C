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

/// \file   writeLUTCollection.C
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Writer for the collection of LUTs for DelphesO2TrackSmearer

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/GeometryContainer.h"

#include <CCDB/BasicCCDBManager.h>

#include <TFile.h>

#include <fairlogger/Logger.h>

void writeLUTCollection(std::string geometryFile = "/home/njacazio/alice/O2Physics/ALICE3/Macros/Configuration/a3geo.ini")
{

  fair::Logger::SetConsoleSeverity(fair::Severity::debug);
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  ccdb.setURL("http://alice-ccdb.cern.ch");
  ccdb.setTimestamp(-1);

  const std::string filename = "ccdb:/Users/j/jekarlss/LookUpTables/NoEloss/el";
  const std::string localFilename = o2::fastsim::GeometryEntry::accessFile(filename, "./.ALICE3/LUTs/", &ccdb, -10);

  o2::delphes::DelphesO2TrackSmearer mSmearer;
  mSmearer.loadTable(11, localFilename.c_str(), true);

  TFile outFile("/tmp/LUTCollection.root", "RECREATE");
  outFile.cd();
  mSmearer.Write();
  outFile.Close();
}
