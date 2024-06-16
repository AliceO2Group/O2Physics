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
/// \file   MetadataHelper.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base of utilities to build advanced tasks
///

#include "Common/Core/MetadataHelper.h"

#include "Framework/InitContext.h"
#include "Framework/RunningWorkflowInfo.h"

MetadataHelper::MetadataHelper()
{
  const std::array<std::string, 5> keyList = {"DataType",
                                              "RecoPassName",
                                              "Run",
                                              "AnchorPassName",
                                              "AnchorProduction"};
  for (const auto& key : keyList) {
    mMetadata[key] = "undefined";
  }
}

void MetadataHelper::initMetadata(o2::framework::ConfigContext const& cfgc)
{
  if (mIsInitialized) {
    LOG(fatal) << "Metadata already initialized. Cannot reinitialize";
  }
  for (const auto& key : mMetadata) {
    const std::string cfgKey = "aod-metadata-" + key.first;
    if (cfgc.options().hasOption(cfgKey.c_str())) {
      mMetadata[key.first] = cfgc.options().get<std::string>(cfgKey.c_str());
      LOG(info) << "Setting metadata " << key.first << " to '" << mMetadata[key.first] << "'";
    }
  }
  mIsInitialized = true;
}

void MetadataHelper::print() const
{
  if (!mIsInitialized) {
    LOG(fatal) << "Metadata not initialized";
  }
  for (const auto& key : mMetadata) {
    LOG(info) << "Metadata " << key.first << ": " << key.second;
  }
}

bool MetadataHelper::isKeyDefined(const std::string& key) const
{
  if (!mIsInitialized) {
    LOG(fatal) << "Metadata not initialized";
  }
  if (mMetadata.find(key) == mMetadata.end()) {
    LOG(fatal) << "Key " << key << " not found in metadata";
  }
  return mMetadata.at(key) != "undefined";
}

bool MetadataHelper::isFullyDefined() const
{
  if (!mIsInitialized) {
    LOG(fatal) << "Metadata not initialized";
  }
  for (const auto& key : mMetadata) {
    if (!isKeyDefined(key.first)) {
      return false;
    }
  }
  return true;
}

std::string MetadataHelper::get(std::string const& key) const
{
  if (!mIsInitialized) {
    LOG(fatal) << "Metadata not initialized";
  }
  if (mMetadata.find(key) == mMetadata.end()) {
    LOG(fatal) << "Key " << key << " not found in metadata";
  }
  return mMetadata.at(key);
}

bool MetadataHelper::isRun3() const
{
  const bool b = (get("Run") == "3");
  LOG(info) << "From metadata this data is from " << (b ? "Run 3" : "Run 2");
  return b;
}

bool MetadataHelper::isMC() const
{
  const bool b = (get("DataType") == "MC");
  LOG(info) << "From metadata this data is from " << (b ? "MC" : "Data");
  return b;
}

bool MetadataHelper::isInitialized() const
{
  if (mIsInitialized) {
    LOG(debug) << "Metadata is initialized";
  } else {
    LOG(debug) << "Metadata is not initialized";
  }
  return mIsInitialized;
}
