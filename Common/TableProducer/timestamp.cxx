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
#include "Common/Core/MetadataHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <RtypesCore.h>

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

struct TimestampTask {
  Produces<aod::Timestamps> timestampTable; /// Table with SOR timestamps produced by the task
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// CCDB manager to access orbit-reset timestamp
  o2::ccdb::CcdbApi ccdb_api;               /// API to access CCDB headers
  Configurable<bool> fatalOnInvalidTimestamp{"fatalOnInvalidTimestamp", false, "Generate fatal error for invalid timestamps"};
  std::map<int, int64_t> mapRunToOrbitReset;                      /// Cache of orbit reset timestamps
  std::map<int, std::pair<int64_t, int64_t>> mapRunToRunDuration; /// Cache of run duration timestamps
  int lastRunNumber = 0;                                          /// Last run number processed
  int64_t orbitResetTimestamp = 0;                                /// Orbit-reset timestamp in us
  std::pair<int64_t, int64_t> runDuration;                        /// Pair of SOR and EOR timestamps

  // Configurables
  Configurable<bool> verbose{"verbose", false, "verbose mode"};
  Configurable<std::string> rct_path{"rct-path", "RCT/Info/RunInformation", "path to the ccdb RCT objects for the SOR timestamps"};
  Configurable<std::string> orbit_reset_path{"orbit-reset-path", "CTP/Calib/OrbitReset", "path to the ccdb orbit-reset objects"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB database"};
  Configurable<int> isRun2MC{"isRun2MC", -1, "Running mode: enable only for Run 2 MC. Timestamps are set to SOR timestamp. Default: -1 (autoset from metadata) 0 (Standard) 1 (Run 2 MC)"};

  void init(o2::framework::InitContext&)
  {
    LOGF(info, "Initializing TimestampTask");
    ccdb->setURL(url.value); // Setting URL of CCDB manager from configuration
    ccdb_api.init(url.value);
    if (!ccdb_api.isHostReachable()) {
      LOGF(fatal, "CCDB host %s is not reacheable, cannot go forward", url.value.data());
    }
    if (isRun2MC.value == -1) {
      if ((!metadataInfo.isRun3()) && metadataInfo.isMC()) {
        isRun2MC.value = 1;
        LOG(info) << "Autosetting the Run2 MC mode based on metadata";
      } else {
        isRun2MC.value = 0;
      }
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
      runDuration = mapRunToRunDuration[runNumber];
    } else { // The run was not requested before: need to acccess CCDB!
      LOGF(debug, "Getting start-of-run and end-of-run timestamps from CCDB");
      runDuration = ccdb->getRunDuration(runNumber, true); /// fatalise if timestamps are not found
      int64_t sorTimestamp = runDuration.first;            // timestamp of the SOR/SOX/STF in ms
      int64_t eorTimestamp = runDuration.second;           // timestamp of the EOR/EOX/ETF in ms

      const bool isUnanchoredRun3MC = runNumber >= 300000 && runNumber < 500000;
      if (isRun2MC.value == 1 || isUnanchoredRun3MC) {
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
        auto ctp = ccdb->getForTimeStamp<std::vector<Long64_t>>(orbit_reset_path.value.data(), eorTimestamp / 2 + sorTimestamp / 2);
        orbitResetTimestamp = (*ctp)[0];
      }

      // Adding the timestamp to the cache map
      std::pair<std::map<int, int64_t>::iterator, bool> check;
      check = mapRunToOrbitReset.insert(std::pair<int, int64_t>(runNumber, orbitResetTimestamp));
      if (!check.second) {
        LOGF(fatal, "Run number %i already existed with a orbit-reset timestamp of %llu", runNumber, check.first->second);
      }
      mapRunToRunDuration[runNumber] = runDuration;
      LOGF(info, "Add new run number %i with orbit-reset timestamp %llu, SOR: %llu, EOR: %llu to cache", runNumber, orbitResetTimestamp, runDuration.first, runDuration.second);
    }

    if (verbose.value) {
      LOGF(info, "Orbit-reset timestamp for run number %i found: %llu us", runNumber, orbitResetTimestamp);
    }
    int64_t timestamp{(orbitResetTimestamp + int64_t(bc.globalBC() * o2::constants::lhc::LHCBunchSpacingNS * 1e-3)) / 1000}; // us -> ms
    if (timestamp < runDuration.first || timestamp > runDuration.second) {
      if (fatalOnInvalidTimestamp.value) {
        LOGF(fatal, "Timestamp %llu us is out of run duration [%llu, %llu] ms", timestamp, runDuration.first, runDuration.second);
      } else {
        LOGF(debug, "Timestamp %llu us is out of run duration [%llu, %llu] ms", timestamp, runDuration.first, runDuration.second);
      }
    }
    timestampTable(timestamp);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{adaptAnalysisTask<TimestampTask>(cfgc)};
}
