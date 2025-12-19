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
/// \file   handleParamBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-02-23
/// \brief  Header file with utilities for handling PID parametrization on CCDB
///

#ifndef COMMON_TOOLS_PID_HANDLEPARAMBASE_H_
#define COMMON_TOOLS_PID_HANDLEPARAMBASE_H_

#include <CCDB/CcdbApi.h>
#include <Framework/Logger.h>

#include <TString.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <map>
#include <string>

// Global executable arguments
namespace bpo = boost::program_options;
bpo::variables_map arguments;             // Command line arguments
o2::ccdb::CcdbApi api;                    // Global CCDB api
unsigned int minRunNumber = 0;            // Starting run validity
unsigned int maxRunNumber = minRunNumber; // Ending run validity
int64_t ccdbTimestamp = 0;                // Timestamp used for the retrieval
int64_t validityStart = 0;                // Initial validity for the object
int64_t validityStop = 0;                 // End validity for the object

std::string timeStampToHReadble(time_t rawtime)
{
  if (rawtime < 0) {
    return std::string(" latest");
  }
  struct tm* dt = new tm();
  char buffer[30];
  rawtime /= 1000;
  localtime_r(&rawtime, dt);
  strftime(buffer, sizeof(buffer), "%H:%M %d-%m %Y", dt);
  return std::string(buffer);
}

// Initializer of the CCDB API
void initCCDBApi()
{
  const auto url = arguments["url"].as<std::string>();
  LOG(info) << "Init CCDB api to URL: " << url;
  api.init(url);
  if (!api.isHostReachable()) {
    LOG(fatal) << "CCDB host " << url << " is not reacheable, cannot go forward";
  }
}

void setStandardOpt(bpo::options_description& options)
{
  options.add_options()(
    "dryrun,D", bpo::value<int>()->default_value(1), "Dryrun mode")(
    "url,u", bpo::value<std::string>()->default_value("http://alice-ccdb.cern.ch"), "URL of the CCDB database e.g. http://ccdb-test.cern.ch:8080 or http://alice-ccdb.cern.ch")(
    "rct-path", bpo::value<std::string>()->default_value("RCT/Info/RunInformation"), "path to the ccdb RCT objects for the SOR/EOR timestamps")(
    "start,s", bpo::value<int64_t>()->default_value(0), "Start timestamp of object validity. If 0 and min-runnumber != 0 it will be set to the run SOR")(
    "stop,S", bpo::value<int64_t>()->default_value(0), "Stop timestamp of object validity. If 0 and max-runnumber != 0 it will be set to the run EOR")(
    "timestamp,T", bpo::value<int64_t>()->default_value(-1), "Timestamp of the object to retrieve, used in alternative to the run number")(
    "min-runnumber,r", bpo::value<unsigned int>()->default_value(0), "Starting run number validity (included) corresponding to the parametrization")(
    "max-runnumber,R", bpo::value<unsigned int>()->default_value(0), "Ending run number validity (included) corresponding to the parametrization. If not specified coincides with min-runnumber")(
    "delete-previous,delete_previous,d", bpo::value<int>()->default_value(0), "Flag to delete previous versions of converter objects in the CCDB before uploading the new one so as to avoid proliferation on CCDB")(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "verbose,v", bpo::value<int>()->default_value(0), "Verbose level 0, 1")(
    "help,h", "Produce help message.");
}

template <typename T>
T* retrieveFromCCDB(const std::string& path,
                    const int64_t timestamp,
                    const std::map<std::string, std::string>& metadata)
{
  std::map<std::string, std::string> headers;
  LOG(info) << "Object " << path << " for timestamp " << timestamp << " -> " << timeStampToHReadble(timestamp);
  headers = api.retrieveHeaders(path, metadata, timestamp);
  LOG(info) << headers.size() << " HEADERS:";
  for (auto const& [key, val] : headers) {
    LOG(info) << "  '" << key << "' -> " << val;
  }

  T* obj = api.retrieveFromTFileAny<T>(path, metadata, timestamp, /*std::map<std::string, std::string>* headers =*/nullptr);
  if (!obj) {
    LOG(fatal) << "Could not get from CCDB for path " << path << " and timestamp " << timestamp << " -> " << timeStampToHReadble(timestamp);
  }
  return obj;
}

template <typename T>
T* retrieveFromCCDB(const std::string& path,
                    const int64_t timestamp)
{
  std::map<std::string, std::string> metadata;
  return retrieveFromCCDB<T>(path, timestamp, metadata);
}

template <typename T>
void storeOnCCDB(const std::string& path,
                 const std::map<std::string, std::string>& metadata,
                 const int64_t& start,
                 const int64_t& stop,
                 const T* obj)
{
  const auto dryrun = arguments["dryrun"].as<int>();
  if (!dryrun) {
    if (arguments["delete-previous"].as<int>()) {
      LOG(info) << "Truncating CCDB path " << path;
      api.truncate(path);
    }
    api.storeAsTFileAny(obj, path, metadata, start, stop);
  } else {
    LOG(info) << "Dryrunning 'api.storeAsTFileAny(obj," << path << ", metadata, " << start << ", " << stop << ")'";
  }
}

void setupTimestamps(int64_t& timestamp,
                     int64_t&,
                     int64_t&)
{
  ccdbTimestamp = arguments["timestamp"].as<int64_t>();
  validityStart = arguments["start"].as<int64_t>();
  validityStop = arguments["stop"].as<int64_t>();
  minRunNumber = arguments["min-runnumber"].as<unsigned int>();
  auto mrun = arguments["max-runnumber"].as<unsigned int>();
  maxRunNumber = mrun > 0 ? mrun : minRunNumber;
  if (minRunNumber > maxRunNumber) {
    LOG(fatal) << "Cannot have `min-runnumber` " << minRunNumber << " > `max-runnumber`" << maxRunNumber;
  }

  auto getSOREOR = [&](const unsigned int runnumber, int64_t& sor, int64_t& eor) {
    std::map<std::string, std::string> metadata, headers;
    const auto rct_path = arguments["rct-path"].as<std::string>();
    const std::string run_path = Form("%s/%i", rct_path.data(), runnumber);

    headers = api.retrieveHeaders(run_path, metadata, -1);

    if (headers.count("STF") == 0) {
      LOGF(warning, "Cannot find STF for run %d in path '%s'. Using SOR as fallback", runnumber, run_path.data());
      if (headers.count("SOR") == 0) {
        LOGF(fatal, "Cannot find SOR in path '%s'.", run_path.data());
      }
      sor = atol(headers["SOR"].c_str());
    } else {
      sor = atol(headers["STF"].c_str());
    }

    if (headers.count("ETF") == 0) {
      LOGF(warning, "Cannot find ETF for run %d in path '%s'. Using EOR as fallback", runnumber, run_path.data());
      if (headers.count("EOR") == 0) {
        LOGF(fatal, "Cannot find EOR in path '%s'.", run_path.data());
      }
      eor = atol(headers["EOR"].c_str());
    } else {
      eor = atol(headers["ETF"].c_str());
    }

    LOG(info) << "Getting timestamp for run " << runnumber << " from CCDB in path " << run_path << " -> STF " << sor << " (" << timeStampToHReadble(sor) << ")"
              << ", ETF " << eor << " (" << timeStampToHReadble(eor) << ")";
  };

  if (minRunNumber != 0) {
    int64_t SOR = 0, EOR = 0;
    getSOREOR(minRunNumber, SOR, EOR);
    timestamp = SOR; // timestamp of the SOR in ms
    LOG(info) << "Setting timestamp of object from run number " << minRunNumber << ": " << timestamp << " -> " << timeStampToHReadble(timestamp);
    if (validityStart == 0) { // Start of validity from first run number
      validityStart = SOR;
      LOG(info) << "Setting validityStart of object from run number " << minRunNumber << ": " << validityStart << " -> " << timeStampToHReadble(validityStart);
      validityStart -= 120000; // add 2 minute margin before start of validity from RCT
      LOG(info) << "Adding 2-minute margin to validityStart: " << validityStart << " -> " << timeStampToHReadble(validityStart);
    }
    if (validityStop == 0) {
      if (minRunNumber != maxRunNumber) {
        getSOREOR(maxRunNumber, SOR, EOR);
        validityStop = EOR;
        LOG(info) << "Setting validityStop of object from run number " << maxRunNumber << ": " << validityStop << " -> " << timeStampToHReadble(validityStop) << " duration " << validityStop - validityStart << " ms, " << (validityStop - validityStart) / 1000 / 3600 << " hours";
      } else {
        validityStop = EOR;
        LOG(info) << "Setting validityStop of object from run number " << minRunNumber << ": " << validityStop << " -> " << timeStampToHReadble(validityStop);
      }
      validityStop += 120000; // add 2 minute margin after end of validity from RCT
      LOG(info) << "Adding 2-minute margin to validityStop: " << validityStop << " -> " << timeStampToHReadble(validityStop);
    }
  }
  if (validityStop == 0) { // Default value for validityStop
    validityStop = 4108971600000;
  }
}

#endif // COMMON_TOOLS_PID_HANDLEPARAMBASE_H_
