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

void MetadataHelper::initMetadata(o2::framework::ConfigContext const& cfgc)
{
  if (cfgc.options().hasOption("aod-metadata-DataType")) {
    aodmetadataDataType = cfgc.options().get<std::string>("aod-metadata-DataType");
    LOG(info) << "Setting metadata DataType to '" << aodmetadataDataType << "'";
  }
  if (cfgc.options().hasOption("aod-metadata-RecoPassName")) {
    aodmetadataRecoPassName = cfgc.options().get<std::string>("aod-metadata-RecoPassName");
    LOG(info) << "Setting metadata RecoPassName to '" << aodmetadataRecoPassName << "'";
  }
  if (cfgc.options().hasOption("aod-metadata-Run")) {
    aodmetadataRun = cfgc.options().get<std::string>("aod-metadata-Run");
    LOG(info) << "Setting metadata Run to '" << aodmetadataRun << "'";
  }
  if (cfgc.options().hasOption("aod-metadata-AnchorPassName")) {
    aodmetadataAnchorPassName = cfgc.options().get<std::string>("aod-metadata-AnchorPassName");
    LOG(info) << "Setting metadata AnchorPassName to '" << aodmetadataAnchorPassName << "'";
  }
  if (cfgc.options().hasOption("aod-metadata-AnchorProduction")) {
    aodmetadataAnchorProduction = cfgc.options().get<std::string>("aod-metadata-AnchorProduction");
    LOG(info) << "Setting metadata AnchorProduction to '" << aodmetadataAnchorProduction << "'";
  }
}

void MetadataHelper::print() const
{
  LOG(info) << "Metadata DataType: " << aodmetadataDataType;
  LOG(info) << "Metadata RecoPassName: " << aodmetadataRecoPassName;
  LOG(info) << "Metadata Run: " << aodmetadataRun;
  LOG(info) << "Metadata AnchorPassName: " << aodmetadataAnchorPassName;
  LOG(info) << "Metadata AnchorProduction: " << aodmetadataAnchorProduction;
}

bool MetadataHelper::isFullyDefined() const
{
  if (aodmetadataDataType == "undefined") {
    return false;
  }
  if (aodmetadataRecoPassName == "undefined") {
    return false;
  }
  if (aodmetadataRun == "undefined") {
    return false;
  }
  if (aodmetadataAnchorPassName == "undefined") {
    return false;
  }
  if (aodmetadataAnchorProduction == "undefined") {
    return false;
  }
  return true;
}

int MetadataHelper::isRun3() const
{
  if (aodmetadataRun == "undefined") {
    LOG(error) << "Metadata Run is undefined";
    return -1;
  }
  const bool b = (aodmetadataRun == "3");
  LOG(info) << "From metadata this data is from " << (b ? "Run 3" : "Run 2");
  return b;
}

int MetadataHelper::isMC() const
{
  if (aodmetadataDataType == "undefined") {
    LOG(error) << "Metadata DataType is undefined";
    return -1;
  }
  const bool b = (aodmetadataDataType == "MC");
  LOG(info) << "From metadata this data is from " << (b ? "MC" : "Data");
  return b;
}

bool MetadataHelper::isInitialized() const
{
  if (mIsInitialized) {
    LOG(info) << "Metadata is initialized";
  } else {
    LOG(info) << "Metadata is not initialized";
  }
  return mIsInitialized;
}
