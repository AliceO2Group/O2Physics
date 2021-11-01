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
/// \file   timestamp.cxx
/// \author Nicolò Jacazio
/// \since  2020-06-22
/// \brief  A task to fill the timestamp table from run number.
///         Uses headers from CCDB
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include <CCDB/BasicCCDBManager.h>
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsRaw/HBFUtils.h"
#include <map>

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

struct TimestampTask {
  Produces<aod::Timestamps> timestampTable;    /// Table with SOR timestamps produced by the task
  Service<o2::ccdb::BasicCCDBManager> ccdb;    /// Object manager in CCDB
  o2::ccdb::CcdbApi ccdb_api;                  /// API to access CCDB
  std::map<int, int>* mapStartOrbit = nullptr; /// Map of the starting orbit for the run
  std::map<int, long> mapRunToTimestamp;       /// Cache of processed run numbers
  int lastRunNumber = 0;                       /// Last run number processed
  long runNumberTimeStamp = 0;                 /// Timestamp of the run number, used in the process function to work out the timestamp of the BC
  uint32_t initialOrbit = 0;                   /// Index of the first orbit of the run number, used in the process function to evaluate the offset with respect to the starting of the run
  static constexpr uint16_t initialBC = 0;     /// Index of the initial bc, exact bc number not relevant due to ms precision of timestamps
  InteractionRecord initialIR;                 /// Initial interaction record, used to compute the delta with respect to the start of the run

  // Configurables
  Configurable<bool> verbose{"verbose", false, "verbose mode"};
  Configurable<std::string> rct_path{"rct-path", "RCT/RunInformation/", "path to the ccdb RCT objects for the SOR timestamps"};
  Configurable<std::string> start_orbit_path{"start-orbit-path", "GRP/StartOrbit", "path to the ccdb SOR orbit objects"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB database"};
  Configurable<bool> isRun2MC{"isRun2MC", false, "Running mode: enable only for Run 2 MC. The timestamp of the BC is computed from initialBC and initialOrbit and does not use the global BC."};

  void init(o2::framework::InitContext&)
  {
    LOGF(info, "Initializing TimestampTask");
    ccdb->setURL(url.value); // Setting URL of CCDB manager from configuration
    LOGF(debug, "Getting SOR orbit map from CCDB url '%s' path '%s'", url.value, start_orbit_path.value);
    mapStartOrbit = ccdb->get<std::map<int, int>>(start_orbit_path.value);
    if (!mapStartOrbit) {
      LOGF(fatal, "Cannot find map of SOR orbits in CCDB in path %s", start_orbit_path.value.data());
    }
    ccdb_api.init(url.value);
    if (!ccdb_api.isHostReachable()) {
      LOGF(fatal, "CCDB host %s is not reacheable, cannot go forward", url.value.data());
    }
  }

  void makeInitialOrbit(aod::BC const& bunchCrossing)
  {
    if (!mapStartOrbit->count(bunchCrossing.runNumber())) {
      LOGF(fatal, "Cannot find run %i in mapStartOrbit map", bunchCrossing.runNumber());
    }
    initialOrbit = mapStartOrbit->at(bunchCrossing.runNumber());
    initialIR.bc = initialBC;
    initialIR.orbit = initialOrbit;
    // Setting lastCall
    LOGF(debug, "Setting the last call of the timestamp for run %i to %llu", bunchCrossing.runNumber(), runNumberTimeStamp);
    lastRunNumber = bunchCrossing.runNumber(); // Setting latest run number information
  }

  void process(aod::BC const& bc)
  {
    // First: we need to set the timestamp from the run number.
    // This is done with caching if the run number of the BC was already processed before
    // If not the timestamp of the run number from BC is queried from CCDB and added to the cache
    if (bc.runNumber() == lastRunNumber) { // The run number coincides to the last run processed
      LOGF(debug, "Using timestamp from last call");
    } else if (mapRunToTimestamp.count(bc.runNumber())) { // The run number was already requested before: getting it from cache!
      LOGF(debug, "Getting timestamp from cache");
      runNumberTimeStamp = mapRunToTimestamp[bc.runNumber()];
      makeInitialOrbit(bc);
    } else { // The run was not requested before: need to acccess CCDB!
      LOGF(debug, "Getting timestamp from CCDB");
      std::map<std::string, std::string> metadata, headers;
      const std::string run_path = Form("%s/%i", rct_path.value.data(), bc.runNumber());
      headers = ccdb_api.retrieveHeaders(run_path, metadata, -1);
      if (headers.count("SOR") == 0) {
        LOGF(fatal, "Cannot find run-number to timestamp in path '%s'.", run_path.data());
      }
      runNumberTimeStamp = atol(headers["SOR"].c_str()); // timestamp of the SOR in ms

      // Adding the timestamp to the cache map
      std::pair<std::map<int, long>::iterator, bool> check;
      check = mapRunToTimestamp.insert(std::pair<int, long>(bc.runNumber(), runNumberTimeStamp));
      if (!check.second) {
        LOGF(fatal, "Run number %i already existed with a timestamp of %llu", bc.runNumber(), check.first->second);
      }
      makeInitialOrbit(bc);
      LOGF(info, "Add new run number %i with timestamp %llu to cache", bc.runNumber(), runNumberTimeStamp);
    }

    if (verbose.value) {
      LOGF(info, "Run-number to timestamp found! %i %llu ms", bc.runNumber(), runNumberTimeStamp);
    }
    const uint16_t currentBC = isRun2MC ? initialBC : (bc.globalBC() % o2::constants::lhc::LHCMaxBunches);
    const uint32_t currentOrbit = isRun2MC ? initialOrbit : (bc.globalBC() / o2::constants::lhc::LHCMaxBunches);
    const InteractionRecord currentIR(currentBC, currentOrbit);
    timestampTable(runNumberTimeStamp + (currentIR - initialIR).bc2ns() * 1e-6);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TimestampTask>(cfgc)};
}
