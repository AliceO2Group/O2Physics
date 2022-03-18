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

#include "CCDB/CcdbApi.h"
#include <boost/program_options.hpp>
#include <FairLogger.h>
#include "TFile.h"

// Global executable arguments
namespace bpo = boost::program_options;
bpo::variables_map arguments;

std::string timeStampToHReadble(time_t rawtime)
{
  if (rawtime < 0) {
    return std::string(" latest");
  }
  struct tm* dt;
  char buffer[30];
  rawtime /= 1000;
  dt = localtime(&rawtime);
  strftime(buffer, sizeof(buffer), "%H:%M %d-%m %Y", dt);
  return std::string(buffer);
}

// Global CCDB api
o2::ccdb::CcdbApi api;

template <typename T>
T* retrieveFromCCDB(const std::string path,
                    const long timestamp)
{
  std::map<std::string, std::string> metadata, headers;
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
void storeOnCCDB(const std::string& path,
                 const std::map<std::string, std::string>& metadata,
                 const long& start,
                 const long& stop,
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

void setupTimestamps(long& timestamp,
                     long& start,
                     long& stop)
{

  const auto runnumber = arguments["runnumber"].as<unsigned int>();
  if (runnumber != 0) {
    std::map<std::string, std::string> metadata, headers;
    const auto rct_path = arguments["rct-path"].as<std::string>();
    const std::string run_path = Form("%s/%i", rct_path.data(), runnumber);

    LOG(info) << "Getting timestamp for run " << runnumber << " from CCDB in path " << run_path;
    headers = api.retrieveHeaders(run_path, metadata, -1);
    if (headers.count("SOR") == 0) {
      LOGF(fatal, "Cannot find run-number to timestamp in path '%s'.", run_path.data());
    }
    if (headers.count("SOR") == 0) {
      LOGF(fatal, "Cannot find run-number SOR in path '%s'.", run_path.data());
    }
    if (headers.count("EOR") == 0) {
      LOGF(fatal, "Cannot find run-number EOR in path '%s'.", run_path.data());
    }
    timestamp = atol(headers["SOR"].c_str()); // timestamp of the SOR in ms
    if (start == 0) {

      start = atol(headers["SOR"].c_str()); // timestamp of the SOR in ms
      LOG(info) << "Setting start of object from run number: " << start << " -> " << timeStampToHReadble(start);
    }
    if (stop == 0) {
      stop = atol(headers["EOR"].c_str()); // timestamp of the EOR in ms
      LOG(info) << "Setting stop of object from run number: " << stop << " -> " << timeStampToHReadble(stop);
    }
  }
  if (stop == 0) { //Default value for stop
    stop = 4108971600000;
  }
}
