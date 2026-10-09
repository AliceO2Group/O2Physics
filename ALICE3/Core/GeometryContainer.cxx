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
/// \file   GeometryContainer.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 geometry handling
/// \since  February 13, 2026
///

#include "GeometryContainer.h"

#include "Common/Core/TableHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>

#include <map>
#include <string>
#include <vector>

namespace o2::fastsim
{

bool GeometryContainer::mCleanLutWhenLoaded = true;
void GeometryContainer::init(o2::framework::InitContext& initContext)
{
  std::vector<std::string> detectorConfiguration;
  const bool foundDetectorConfiguration = common::core::getTaskOptionValue(initContext, "on-the-fly-detector-geometry-provider", "detectorConfiguration", detectorConfiguration, false);
  if (!foundDetectorConfiguration) {
    LOG(fatal) << "Could not retrieve detector configuration from OnTheFlyDetectorGeometryProvider task.";
    return;
  }
  LOG(info) << "Size of detector configuration: " << detectorConfiguration.size();

  bool cleanLutWhenLoaded;
  const bool foundCleanLutWhenLoaded = common::core::getTaskOptionValue(initContext, "on-the-fly-detector-geometry-provider", "cleanLutWhenLoaded", cleanLutWhenLoaded, false);
  if (!foundCleanLutWhenLoaded) {
    LOG(fatal) << "Could not retrieve foundCleanLutWhenLoaded option from OnTheFlyDetectorGeometryProvider task.";
    return;
  }
  setLutCleanupSetting(cleanLutWhenLoaded);

  for (std::string& configFile : detectorConfiguration) {
    LOG(info) << "Detector geometry configuration file used: " << configFile;
    addEntry(configFile);
  }
}

void GeometryContainer::addEntry(const std::string& filename)
{
  if (!mCcdb) {
    LOG(fatal) << " --- ccdb is not set";
  }
  mEntries.emplace_back(filename, mCcdb);
}

std::map<std::string, std::string> GeometryEntry::getConfiguration(const std::string& layerName) const
{
  auto it = mConfigurations.find(layerName);
  if (it != mConfigurations.end()) {
    return it->second;
  } else {
    LOG(fatal) << "Layer " << layerName << " not found in geometry configurations.";
    return {};
  }
}

bool GeometryEntry::hasValue(const std::string& layerName, const std::string& key) const
{
  auto layerIt = mConfigurations.find(layerName);
  if (layerIt != mConfigurations.end()) {
    auto keyIt = layerIt->second.find(key);
    return keyIt != layerIt->second.end();
  }
  return false;
}

std::string GeometryEntry::getValue(const std::string& layerName, const std::string& key, bool require) const
{
  auto layer = getConfiguration(layerName);
  auto entry = layer.find(key);
  if (entry != layer.end()) {
    return layer.at(key);
  } else if (require) {
    LOG(fatal) << "Key " << key << " not found in layer " << layerName << " configurations.";
    return "";
  } else {
    return "";
  }
}

void GeometryEntry::replaceValue(const std::string& layerName, const std::string& key, const std::string& value)
{
  if (!hasValue(layerName, key)) { // check that the key exists
    LOG(fatal) << "Key " << key << " does not exist in layer " << layerName << ". Cannot replace value.";
  }
  setValue(layerName, key, value);
}

} // namespace o2::fastsim
