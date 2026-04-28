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
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <string>
#include <vector>

struct OnTheFlyDetectorGeometryProvider {
  o2::framework::HistogramRegistry histos{"Histos", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  o2::framework::Configurable<bool> cleanLutWhenLoaded{"cleanLutWhenLoaded", true, "clean LUTs after being loaded to save disk space"};
  o2::framework::Configurable<std::vector<std::string>> detectorConfiguration{"detectorConfiguration",
                                                                              std::vector<std::string>{"$O2PHYSICS_ROOT/share/alice3/a3geometry_v3.ini"},
                                                                              "Paths of the detector geometry configuration files"};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);
    o2::fastsim::GeometryContainer geometryContainer; // Checking that the geometry files can be accessed and loaded
    geometryContainer.setCcdbManager(ccdb.operator->());
    LOG(info) << "On-the-fly detector geometry provider running.";
    if (detectorConfiguration.value.empty()) {
      LOG(fatal) << "No detector configuration files provided.";
      return;
    }
    int idx = 0;
    for (std::string& configFile : detectorConfiguration.value) {
      LOG(info) << "Loading detector geometry from configuration file: " << configFile;
      histos.add<TH1>(Form("GeometryConfigFile_%d", idx++), configFile.c_str(), o2::framework::HistType::kTH1D, {{1, 0, 1}})->Fill(0.5);
      geometryContainer.addEntry(configFile);
    }

    // First we check that the magnetic field is consistent
    const int nGeometries = geometryContainer.getNumberOfConfigurations();
    const float mMagneticField = geometryContainer.getFloatValue(0, "global", "magneticfield");
    for (int icfg = 0; icfg < nGeometries; ++icfg) {
      const float cfgBfield = geometryContainer.getFloatValue(icfg, "global", "magneticfield");
      if (std::abs(cfgBfield - mMagneticField) > 1e-3) {
        LOG(fatal) << "Inconsistent magnetic field values between configurations 0 and " << icfg << ": " << mMagneticField << " vs " << cfgBfield;
      }
    }
    LOG(info) << "Initialization completed";
  }

  void process(o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles)
  {
    LOG(debug) << "On-the-fly detector geometry provider processing " << mcCollisions.size() << " collisions and " << mcParticles.size() << " particles.";
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
