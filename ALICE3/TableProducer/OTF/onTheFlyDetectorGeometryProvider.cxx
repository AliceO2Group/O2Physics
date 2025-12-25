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

/// \file onTheFlyDetectorGeometryProvider.cxx
///
/// \brief Reader and and provider of the detector geometry for on-the-fly simulation
///
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, Universita del Piemonte Orientale (IT)
///

#include "ALICE3/Core/FastTracker.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <string>
#include <vector>

struct OnTheFlyDetectorGeometryProvider {

  o2::framework::Configurable<std::vector<std::string>> detectorConfiguration{"detectorConfiguration",
                                                                              std::vector<std::string>{"$O2PHYSICS_ROOT/share/alice3/a3geometry_v3.ini"},
                                                                              "Paths of the detector geometry configuration files"};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);
  }
  void run(o2::framework::ProcessingContext& pc)
  {
    o2::fastsim::GeometryContainer geometryContainer; // Checking that the geometry files can be accessed and loaded
    LOG(info) << "On-the-fly detector geometry provider running.";
    if (detectorConfiguration.value.empty()) {
      LOG(fatal) << "No detector configuration files provided.";
      return;
    }
    int idx = 0;
    for (auto& configFile : detectorConfiguration.value) {
      LOG(info) << "Loading detector geometry from configuration file: " << configFile;
      // If the filename starts with ccdb: then take the file from the ccdb
      if (configFile.rfind("ccdb:", 0) == 0) {
        std::string ccdbPath = configFile.substr(5); // remove "ccdb:" prefix
        const std::string outPath = "/tmp/DetGeo/";
        configFile = Form("%s/%s/snapshot.root", outPath.c_str(), ccdbPath.c_str());
        std::ifstream checkFile(configFile); // Check if file already exists
        if (!checkFile.is_open()) {          // File does not exist, retrieve from CCDB
          LOG(info) << " --- CCDB source detected for detector geometry " << configFile;
          std::map<std::string, std::string> metadata;
          ccdb->getCCDBAccessor().retrieveBlob(ccdbPath, outPath, metadata, 1);
          LOG(info) << " --- Now retrieving geometry configuration from CCDB to: " << configFile;
        } else { // File exists, proceed to load
          LOG(info) << " --- Geometry configuration file already exists: " << configFile << ". Skipping download.";
          checkFile.close();
        }
        detectorConfiguration.value[idx] = configFile; // Update the filename to the local file
      }
      geometryContainer.addEntry(configFile);
      idx++;
    }
    pc.services().get<o2::framework::ControlService>().endOfStream();
    pc.services().get<o2::framework::ControlService>().readyToQuit(o2::framework::QuitRequest::Me);
  }
};

// #define VERIFY_ALICE3
#ifdef VERIFY_ALICE3
struct OnTheFlyDetectorGeometryUser {
  void init(o2::framework::InitContext& initContext)
  {
    o2::fastsim::GeometryContainer geometryContainer; // Checking that the configuration can be inherited
    geometryContainer.init(initContext);
  }
  void run(o2::framework::ProcessingContext& pc)
  {
    pc.services().get<o2::framework::ControlService>().endOfStream();
    pc.services().get<o2::framework::ControlService>().readyToQuit(o2::framework::QuitRequest::Me);
  }
};
#endif

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec spec;
  spec.push_back(adaptAnalysisTask<OnTheFlyDetectorGeometryProvider>(cfgc));
#ifdef VERIFY_ALICE3
  spec.push_back(adaptAnalysisTask<OnTheFlyDetectorGeometryUser>(cfgc));
#endif
  return spec;
}
