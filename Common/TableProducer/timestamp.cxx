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
#include <vector>
#include <map>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsRaw/HBFUtils.h"

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

struct TimestampTask {
  Produces<aod::Timestamps> timestampTable;  /// Table with SOR timestamps produced by the task
  Service<o2::ccdb::BasicCCDBManager> ccdb;  /// CCDB manager to access orbit-reset timestamp
  o2::ccdb::CcdbApi ccdb_api;                /// API to access CCDB headers
  std::map<int, int64_t> mapRunToOrbitReset; /// Cache of orbit reset timestamps
  int lastRunNumber = 0;                     /// Last run number processed
  int64_t orbitResetTimestamp = 0;           /// Orbit-reset timestamp in us

  // Configurables
  Configurable<bool> verbose{"verbose", false, "verbose mode"};
  Configurable<std::string> rct_path{"rct-path", "RCT/Info/RunInformation", "path to the ccdb RCT objects for the SOR timestamps"};
  Configurable<std::string> orbit_reset_path{"orbit-reset-path", "CTP/Calib/OrbitReset", "path to the ccdb orbit-reset objects"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB database"};
  Configurable<bool> isRun2MC{"isRun2MC", false, "Running mode: enable only for Run 2 MC. Timestamps are set to SOR timestamp"};

  void init(o2::framework::InitContext&)
  {
    LOGF(info, "Initializing TimestampTask");
    ccdb->setURL(url.value); // Setting URL of CCDB manager from configuration
    ccdb_api.init(url.value);
    if (!ccdb_api.isHostReachable()) {
      LOGF(fatal, "CCDB host %s is not reacheable, cannot go forward", url.value.data());
    }
  }

  void process(aod::BC const& bc)
  {
    int runNumber = bc.runNumber();
    // We need to set the orbit-reset timestamp for the run number.
    // This is done with caching if the run number was already processed before.
    // If not the orbit-reset timestamp for the run number is queried from CCDB and added to the cache
    if (runNumber == lastRunNumber) { // The run number coincides to the last run processed
      LOGF(debug, "Using orbit-reset timestamp from last call");
    } else if (mapRunToOrbitReset.count(runNumber)) { // The run number was already requested before: getting it from cache!
      LOGF(debug, "Getting orbit-reset timestamp from cache");
      orbitResetTimestamp = mapRunToOrbitReset[runNumber];
    } else { // The run was not requested before: need to acccess CCDB!
      LOGF(debug, "Getting start-of-run and end-of-run timestamps from CCDB");
      auto timestamps = ccdb->getRunDuration(runNumber, true); /// fatalise if timestamps are not found
      int64_t sorTimestamp = timestamps.first;                 // timestamp of the SOR in ms
      int64_t eorTimestamp = timestamps.second;                // timestamp of the EOR in ms

      bool isUnanchoredRun3MC = runNumber >= 300000 && runNumber < 500000;
      if (isRun2MC || isUnanchoredRun3MC) {
        // isRun2MC: bc/orbit distributions are not simulated in Run2 MC. All bcs are set to 0.
        // isUnanchoredRun3MC: assuming orbit-reset is done in the beginning of each run
        // Setting orbit-reset timestamp to start-of-run timestamp
        orbitResetTimestamp = sorTimestamp * 1000; // from ms to us
      } else if (runNumber < 300000) {             // Run 2
        LOGF(debug, "Getting orbit-reset timestamp using start-of-run timestamp from CCDB");
        auto ctp = ccdb->getForTimeStamp<std::vector<Long64_t>>(orbit_reset_path.value.data(), sorTimestamp);
        orbitResetTimestamp = (*ctp)[0];
      } else {
        // sometimes orbit is reset after SOR. Using EOR timestamps for orbitReset query is more reliable
        LOGF(debug, "Getting orbit-reset timestamp using end-of-run timestamp from CCDB");
        auto ctp = ccdb->getForTimeStamp<std::vector<Long64_t>>(orbit_reset_path.value.data(), eorTimestamp);
        orbitResetTimestamp = (*ctp)[0];
      }

      // Adding the timestamp to the cache map
      std::pair<std::map<int, int64_t>::iterator, bool> check;
      check = mapRunToOrbitReset.insert(std::pair<int, int64_t>(runNumber, orbitResetTimestamp));
      if (!check.second) {
        LOGF(fatal, "Run number %i already existed with a orbit-reset timestamp of %llu", runNumber, check.first->second);
      }
      LOGF(info, "Add new run number %i with orbit-reset timestamp %llu to cache", runNumber, orbitResetTimestamp);
    }

    if (verbose.value) {
      LOGF(info, "Orbit-reset timestamp for run number %i found: %llu us", runNumber, orbitResetTimestamp);
    }

    timestampTable((orbitResetTimestamp + int64_t(bc.globalBC() * o2::constants::lhc::LHCBunchSpacingNS * 1e-3)) / 1000); // us -> ms
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TimestampTask>(cfgc)};
}
