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
/// \file   GeometryContainer.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 geometry handling
/// \since  February 13, 2026
///

#ifndef ALICE3_CORE_GEOMETRYCONTAINER_H_
#define ALICE3_CORE_GEOMETRYCONTAINER_H_

#include "ALICE3/Core/ConfigurationParser.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>

#include <map>
#include <string>
#include <vector>

namespace o2::fastsim
{

// A container for the geometry info
struct GeometryEntry {
  // Default constructor
  GeometryEntry() = default;
  explicit GeometryEntry(const std::string& filename, o2::ccdb::BasicCCDBManager* ccdb = nullptr)
  {
    mFileName = ConfigurationParser::accessFile(filename, "./.ALICE3/Configuration/", ccdb);
    mConfigurations = ConfigurationParser::parseTEnvConfiguration(mFileName, mLayerNames);
    LOG(info) << "Loaded geometry configuration from file: " << mFileName << " with " << mLayerNames.size() << " layers.";
    if (mLayerNames.empty()) {
      LOG(warning) << "No layers found in geometry configuration file: " << filename;
    }
  }

  std::map<std::string, std::map<std::string, std::string>> getConfigurations() const { return mConfigurations; }
  std::map<std::string, std::string> getConfiguration(const std::string& layerName) const;
  std::vector<std::string> getLayerNames() const { return mLayerNames; }
  bool hasValue(const std::string& layerName, const std::string& key) const;
  std::string getValue(const std::string& layerName, const std::string& key, bool require = true) const;
  void setValue(const std::string& layerName, const std::string& key, const std::string& value) { mConfigurations[layerName][key] = value; }
  void replaceValue(const std::string& layerName, const std::string& key, const std::string& value);
  float getFloatValue(const std::string& layerName, const std::string& key) const { return std::stof(getValue(layerName, key)); }
  int getIntValue(const std::string& layerName, const std::string& key) const { return std::stoi(getValue(layerName, key)); }

 private:
  std::string mFileName;                                                     // Filename of the geometry
  std::map<std::string, std::map<std::string, std::string>> mConfigurations; // Layer configurations
  std::vector<std::string> mLayerNames;                                      // Ordered names of the layers
};

class GeometryContainer
{
 public:
  GeometryContainer() = default;
  virtual ~GeometryContainer() = default;

  /**
   * @brief Initializes the GeometryContainer by retrieving the list of geometry configuration files from the OnTheFlyDetectorGeometryProvider task options and parsing them to fill the container entries.
   **/
  void init(o2::framework::InitContext& initContext);

  // Add a geometry entry from a configuration file
  void addEntry(const std::string& filename);
  static void setLutCleanupSetting(const bool cleanLutWhenLoaded) { mCleanLutWhenLoaded = cleanLutWhenLoaded; }
  void setCcdbManager(o2::ccdb::BasicCCDBManager* mgr) { mCcdb = mgr; }

  // Getters
  int getNumberOfConfigurations() const { return mEntries.size(); }
  const std::vector<GeometryEntry>& getEntries() const { return mEntries; }
  const GeometryEntry& getEntry(const int id) const { return mEntries.at(id); }
  GeometryEntry getGeometryEntry(const int id) const { return mEntries.at(id); }
  static bool cleanLutWhenLoaded() { return mCleanLutWhenLoaded; }

  // Get configuration maps
  std::map<std::string, std::map<std::string, std::string>> getConfigurations(const int id) const { return mEntries.at(id).getConfigurations(); }
  std::map<std::string, std::string> getConfiguration(const int id, const std::string& layerName) const { return mEntries.at(id).getConfiguration(layerName); }

  // Get specific values
  std::string getValue(const int id, const std::string& layerName, const std::string& key, bool require = true) const { return mEntries.at(id).getValue(layerName, key, require); }
  float getFloatValue(const int id, const std::string& layerName, const std::string& key) const { return mEntries.at(id).getFloatValue(layerName, key); }

 private:
  static bool mCleanLutWhenLoaded; // Whether to clean the LUT when loading a new geometry configuration
  std::vector<GeometryEntry> mEntries;
  o2::ccdb::BasicCCDBManager* mCcdb = nullptr;
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_GEOMETRYCONTAINER_H_
