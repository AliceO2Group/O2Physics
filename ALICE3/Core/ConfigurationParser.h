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
/// \file   ConfigurationParser.h
/// \brief  Utilities to access and parse TEnv configuration files for the ALICE3 fast simulation
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
///

#ifndef ALICE3_CORE_CONFIGURATIONPARSER_H_
#define ALICE3_CORE_CONFIGURATIONPARSER_H_

#include <CCDB/BasicCCDBManager.h>

#include <map>
#include <string>
#include <vector>

namespace o2::fastsim
{

class ConfigurationParser
{
 public:
  /**
   * @brief Parses a TEnv configuration file with keys of the form "<entry>.<parameter>" and returns the key-value pairs split per entry
   * @param filename Path to the TEnv configuration file
   * @param entries Vector to store the order of the entries as they appear in the file
   * @return A map where each key is an entry name and the value is another map of key-value pairs for that entry
   */
  static std::map<std::string, std::map<std::string, std::string>> parseTEnvConfiguration(std::string& filename, std::vector<std::string>& entries);

  /**
   * @brief Accesses a file given its path, which can be either a local path or a ccdb path (starting with "ccdb:"). In the first case it returns the local path, in the second it retrieves the file from ccdb and returns the local path to the retrieved file.
   * @param path The path to the file, either local or ccdb (starting with "ccdb:")
   * @param downloadPath The local path where to download the file if it's a ccdb path. Default is "/tmp/GeometryContainer/"
   * @param ccdb Pointer to the CCDB manager to use for retrieving the file if it's a ccdb path. Must be set when path is a ccdb path.
   * @param timeoutSeconds If positive, then this function will wait for these seconds after download before removing the downloaded file.
   * @return The local path to the file, either the original local path or the path to the retrieved file from ccdb
   */
  static std::string accessFile(const std::string& path, const std::string& downloadPath = "/tmp/GeometryContainer/", o2::ccdb::BasicCCDBManager* ccdb = nullptr, int timeoutSeconds = 0);
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_CONFIGURATIONPARSER_H_
