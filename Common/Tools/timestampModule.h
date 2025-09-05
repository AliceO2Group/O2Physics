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

#ifndef COMMON_TOOLS_TIMESTAMPMODULE_H_
#define COMMON_TOOLS_TIMESTAMPMODULE_H_

#include <CommonConstants/LHCConstants.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace o2
{
namespace common
{
namespace timestamp
{

// timestamp configurables
struct timestampConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "timestamp";
  o2::framework::Configurable<bool> verbose{"verbose", false, "verbose mode"};
  o2::framework::Configurable<bool> fatalOnInvalidTimestamp{"fatalOnInvalidTimestamp", false, "Generate fatal error for invalid timestamps"};
  o2::framework::Configurable<std::string> rct_path{"rct-path", "RCT/Info/RunInformation", "path to the ccdb RCT objects for the SOR timestamps"};
  o2::framework::Configurable<std::string> orbit_reset_path{"orbit-reset-path", "CTP/Calib/OrbitReset", "path to the ccdb orbit-reset objects"};
  o2::framework::Configurable<int> isRun2MC{"isRun2MC", -1, "Running mode: enable only for Run 2 MC. Timestamps are set to SOR timestamp. Default: -1 (autoset from metadata) 0 (Standard) 1 (Run 2 MC)"}; // o2-linter: disable=name/configurable (temporary fix)
};

//__________________________________________
// time stamp module
//
// class to acquire time stamps to be used in
// modular (plugin) fashion

class TimestampModule
{
 public:
  TimestampModule()
  {
    // constructor: initialize at defaults
    lastRunNumber = 0;
    orbitResetTimestamp = 0;
  };

  o2::common::timestamp::timestampConfigurables timestampOpts;

  // objects necessary during processing
  std::map<int, int64_t> mapRunToOrbitReset;                      /// Cache of orbit reset timestamps
  std::map<int, std::pair<int64_t, int64_t>> mapRunToRunDuration; /// Cache of run duration timestamps
  int lastRunNumber;                                              /// Last run number processed
  int64_t orbitResetTimestamp;                                    /// Orbit-reset timestamp in us
  std::pair<int64_t, int64_t> runDuration;                        /// Pair of SOR and EOR timestamps

  template <typename TTimestampOpts, typename TMetadatahelper>
  void init(TTimestampOpts const& external_timestampOpts, TMetadatahelper const& metadataInfo)
  {
    timestampOpts = external_timestampOpts;

    if (timestampOpts.isRun2MC.value == -1) {
      if ((!metadataInfo.isRun3()) && metadataInfo.isMC()) {
        timestampOpts.isRun2MC.value = 1;
        LOG(info) << "Autosetting the Run2 MC mode based on metadata";
      } else {
        timestampOpts.isRun2MC.value = 0;
      }
    }
  }

  template <typename TBCs, typename Tccdb, typename TTimestampBuffer, typename TCursor>
  void process(TBCs const& bcs, Tccdb const& ccdb, TTimestampBuffer& timestampbuffer, TCursor& timestampTable)
  {
    timestampbuffer.clear();
    for (auto const& bc : bcs) {
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

        // clear cache to prevent interference with orbit reset queries from other code
        // FIXME this should not have been a problem, to be investigated
        ccdb->clearCache(timestampOpts.orbit_reset_path.value.data());

        const bool isUnanchoredRun3MC = runNumber >= 300000 && runNumber < 500000;
        if (timestampOpts.isRun2MC.value == 1 || isUnanchoredRun3MC) {
          // isRun2MC: bc/orbit distributions are not simulated in Run2 MC. All bcs are set to 0.
          // isUnanchoredRun3MC: assuming orbit-reset is done in the beginning of each run
          // Setting orbit-reset timestamp to start-of-run timestamp
          orbitResetTimestamp = sorTimestamp * 1000; // from ms to us
        } else if (runNumber < 300000) {             // Run 2
          LOGF(debug, "Getting orbit-reset timestamp using start-of-run timestamp from CCDB");
          auto ctp = ccdb->template getSpecific<std::vector<int64_t>>(timestampOpts.orbit_reset_path.value.data(), sorTimestamp);
          orbitResetTimestamp = (*ctp)[0];
        } else {
          // sometimes orbit is reset after SOR. Using EOR timestamps for orbitReset query is more reliable
          LOGF(debug, "Getting orbit-reset timestamp using end-of-run timestamp from CCDB");
          auto ctp = ccdb->template getSpecific<std::vector<int64_t>>(timestampOpts.orbit_reset_path.value.data(), eorTimestamp / 2 + sorTimestamp / 2);
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

      if (timestampOpts.verbose) {
        LOGF(info, "Orbit-reset timestamp for run number %i found: %llu us", runNumber, orbitResetTimestamp);
      }
      int64_t timestamp{(orbitResetTimestamp + int64_t(bc.globalBC() * o2::constants::lhc::LHCBunchSpacingNS * 1e-3)) / 1000}; // us -> ms
      if (timestamp < runDuration.first || timestamp > runDuration.second) {
        if (timestampOpts.fatalOnInvalidTimestamp.value) {
          LOGF(fatal, "Timestamp %llu us is out of run duration [%llu, %llu] ms", timestamp, runDuration.first, runDuration.second);
        } else {
          LOGF(debug, "Timestamp %llu us is out of run duration [%llu, %llu] ms", timestamp, runDuration.first, runDuration.second);
        }
      }
      timestampbuffer.push_back(timestamp); // for buffering purposes
      timestampTable(timestamp);
    }
  }
};

} // namespace timestamp
} // namespace common
} // namespace o2

#endif // COMMON_TOOLS_TIMESTAMPMODULE_H_
