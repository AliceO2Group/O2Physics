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
/// \file   TableHelper.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base of utilities to build advanced tasks for the metadata
///

#ifndef COMMON_CORE_METADATAHELPER_H_
#define COMMON_CORE_METADATAHELPER_H_

#include <string>

#include "Framework/ConfigContext.h"

struct MetadataHelper {
  std::string aodmetadataDataType = "undefined";
  std::string aodmetadataRecoPassName = "undefined";
  std::string aodmetadataRun = "undefined";
  std::string aodmetadataAnchorPassName = "undefined";
  std::string aodmetadataAnchorProduction = "undefined";

  // Function to initialize the metadata from the configuration context
  void initMetadata(o2::framework::ConfigContext const& cfgc);
  void print() const;

  bool isFullyDefined() const;

  /// Function to check if the data is from Run 3
  /// @return 1 if the run is a run 3 run, 0 if it is not, -1 if it is not defined
  int isRun3() const;

  /// Function to check if the data is from MC
  /// @return 1 if the data is from MC, 0 if it is not, -1 if it is not defined
  int isMC() const;

  /// @brief Function to check if the data has been correctly initialized
  /// @return true if the data has been initialized, false otherwise
  bool isInitialized() const;

 private:
  bool mIsInitialized = false;
};

#endif // COMMON_CORE_METADATAHELPER_H_
