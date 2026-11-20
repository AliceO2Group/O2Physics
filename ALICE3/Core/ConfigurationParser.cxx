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
/// \file   ConfigurationParser.cxx
/// \brief  Utilities to access and parse TEnv configuration files for the ALICE3 fast simulation
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
///

#include "ConfigurationParser.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/Logger.h>

#include <TEnv.h>
#include <THashList.h>
#include <TString.h>
#include <TSystem.h>

#include <sys/file.h>
#include <sys/stat.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace o2::fastsim
{

std::map<std::string, std::map<std::string, std::string>> ConfigurationParser::parseTEnvConfiguration(std::string& filename, std::vector<std::string>& entries)
{
  std::map<std::string, std::map<std::string, std::string>> configMap;
  filename = gSystem->ExpandPathName(filename.c_str());
  LOG(info) << "Parsing TEnv configuration file: " << filename;
  TEnv env(filename.c_str());
  THashList* table = env.GetTable();
  entries.clear();
  for (int i = 0; i < table->GetEntries(); ++i) {
    const std::string key = table->At(i)->GetName();
    // key should contain exactly one dot
    if (key.find('.') == std::string::npos || key.find('.') != key.rfind('.')) {
      LOG(fatal) << "Key " << key << " does not contain exactly one dot";
      continue;
    }
    const std::string firstPart = key.substr(0, key.find('.'));
    if (std::find(entries.begin(), entries.end(), firstPart) == entries.end()) {
      entries.push_back(firstPart);
    }
  }
  env.Print();
  // Entries
  for (const auto& entry : entries) {
    LOG(info) << " Reading entry " << entry;
    for (int i = 0; i < table->GetEntries(); ++i) {
      const std::string key = table->At(i)->GetName();
      if (key.find(entry + ".") == 0) {
        const std::string paramName = key.substr(key.find('.') + 1);
        const std::string value = env.GetValue(key.c_str(), "");
        configMap[entry][paramName] = value;
      }
    }
  }
  return configMap;
}

std::string ConfigurationParser::accessFile(const std::string& path, const std::string& downloadPath, o2::ccdb::BasicCCDBManager* ccdb, int timeoutSeconds)
{
  if (path.rfind("ccdb:", 0) == 0) {
    const std::string ccdbPath = path.substr(5); // remove "ccdb:" prefix
    const std::string localPath = Form("%s/%s/snapshot.root", downloadPath.c_str(), ccdbPath.c_str());
    const std::string lockFile = localPath + ".lock";
    const std::string doneFile = localPath + ".done";

    // Create directory structure if it doesn't exist
    std::string dirPath = localPath.substr(0, localPath.find_last_of('/'));
    gSystem->mkdir(dirPath.c_str(), true);

    // Check if file is already fully downloaded
    struct stat buffer;
    if (stat(doneFile.c_str(), &buffer) == 0) {
      LOG(info) << " --- Geometry configuration file already exists: " << localPath << ". Skipping download.";
      return localPath;
    }

    // Acquire file lock for inter-process synchronization
    int lockFd = open(lockFile.c_str(), O_CREAT | O_RDWR, 0666);
    if (lockFd == -1) {
      LOG(error) << " --- Failed to create lock file: " << lockFile;
      return localPath;
    }

    // Try to acquire exclusive lock (non-blocking)
    LOG(info) << " --- Attempting to acquire lock for: " << localPath;
    int lockResult = flock(lockFd, LOCK_EX | LOCK_NB);

    if (lockResult == -1 && errno == EWOULDBLOCK) {
      // Lock is held by another process - wait up to 10 minutes for download to complete
      LOG(info) << " --- Lock is held by another process. Waiting for download to complete (up to 10 minutes)...";
      close(lockFd);

      const auto startTime = std::chrono::steady_clock::now();
      const auto timeout = std::chrono::minutes(10);
      const auto checkInterval = std::chrono::seconds(5);

      while (true) {
        // Check if download is complete
        if (stat(doneFile.c_str(), &buffer) == 0) {
          LOG(info) << " --- Geometry configuration file was downloaded by another process: " << localPath;
          return localPath;
        }

        // Check timeout
        auto elapsed = std::chrono::steady_clock::now() - startTime;
        if (elapsed >= timeout) {
          LOG(fatal) << " --- Timeout waiting for file download: " << localPath << ". Waited for 10 minutes.";
          return localPath;
        }

        // Wait before checking again
        std::this_thread::sleep_for(checkInterval);
      }
    } else if (lockResult == -1) {
      LOG(error) << " --- Failed to acquire lock for: " << lockFile;
      close(lockFd);
      return localPath;
    }

    // Lock acquired successfully - double-check if file was downloaded while we were trying
    if (stat(doneFile.c_str(), &buffer) == 0) {
      LOG(info) << " --- Geometry configuration file was downloaded by another process: " << localPath;
      flock(lockFd, LOCK_UN);
      close(lockFd);
      return localPath;
    }

    // File does not exist, retrieve from CCDB
    LOG(info) << " --- CCDB source detected for " << path;
    std::map<std::string, std::string> metadata;
    bool status = ccdb->getCCDBAccessor().retrieveBlob(ccdbPath, downloadPath, metadata, 1);
    if (!status) {
      flock(lockFd, LOCK_UN);
      close(lockFd);
      LOG(fatal) << " --- Failed to retrieve configuration from CCDB for path: " << ccdbPath;
      return "";
    }
    LOG(info) << " --- Retrieved configuration from CCDB to: " << localPath;

    // Verify the integrity of the downloaded file
    if (stat(localPath.c_str(), &buffer) != 0) {
      flock(lockFd, LOCK_UN);
      close(lockFd);
      LOG(fatal) << " --- Downloaded file does not exist or is not accessible: " << localPath;
      return "";
    }
    if (buffer.st_size == 0) {
      flock(lockFd, LOCK_UN);
      close(lockFd);
      LOG(fatal) << " --- Downloaded file is empty: " << localPath;
      return "";
    }
    LOG(info) << " --- File integrity verified. Size: " << buffer.st_size << " bytes";

    // Create done marker file to indicate successful download
    std::ofstream doneMarker(doneFile);
    doneMarker.close();

    // Release lock
    flock(lockFd, LOCK_UN);
    close(lockFd);

    // If timeout is specified, schedule file deletion after timeout
    if (timeoutSeconds > 0) {
      LOG(info) << " --- Deleting configuration file after timeout: " << localPath;
      std::thread deletionThread([localPath, doneFile, timeoutSeconds]() {
        LOG(info) << " --- Operating deletion of configuration file after timeout: " << localPath;
        std::this_thread::sleep_for(std::chrono::seconds(timeoutSeconds));
        if (std::remove(localPath.c_str()) == 0) {
          LOG(info) << " --- File deleted successfully: " << localPath;
        } else {
          LOG(warning) << " --- Failed to delete file: " << localPath;
        }
        // Also remove the done marker file
        if (std::remove(doneFile.c_str()) == 0) {
          LOG(info) << " --- Done marker deleted: " << doneFile;
        }
      });
      deletionThread.detach();
    }

    return localPath;
  }
  return path;
}

} // namespace o2::fastsim
