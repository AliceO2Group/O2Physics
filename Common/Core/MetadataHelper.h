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
/// \brief  Utility to handle the metadata from the AOD
///

#ifndef COMMON_CORE_METADATAHELPER_H_
#define COMMON_CORE_METADATAHELPER_H_

#include <string>
#include <map>
#include "Framework/ConfigContext.h"

struct MetadataHelper {
  /// @brief Constructor for the MetadataHelper. Defines the all the metadata keys that will be looked for and accessible
  MetadataHelper();

  /// @brief Function to initialize the metadata from the configuration context.
  /// @param cfgc the configuration context
  void initMetadata(o2::framework::ConfigContext const& cfgc);

  /// @brief Function to print the metadata
  void print() const;

  /// Function to check if the metadata is fully defined
  /// @return true if the metadata is fully defined, false otherwise
  bool isFullyDefined() const;

  /// Function to check if the data is from Run 3
  /// @return 1 if the run is a run 3 run, 0 if it is not, -1 if it is not defined
  bool isRun3() const;

  /// Function to check if the data is from MC
  /// @return 1 if the data is from MC, 0 if it is not, -1 if it is not defined
  bool isMC() const;

  /// @brief Function to check if the data has been correctly initialized
  /// @return true if the data has been initialized, false otherwise
  bool isInitialized() const;

  /// @brief Function to get the metadata value for a given key
  /// @param key the key of the metadata
  /// @return the value of the metadata. Throws an exception if the key is not found
  std::string get(const std::string& key) const;

  /// @brief  Function to check if a key is defined in the metadata
  /// @param key the key to check
  /// @return true if the key is defined, false otherwise. Throws an exception if the key is not found
  bool isKeyDefined(const std::string& key) const;

 private:
  std::map<std::string, std::string> mMetadata; /// < The metadata map
  bool mIsInitialized = false;                  /// < Flag to check if the metadata has been initialized
};

#endif // COMMON_CORE_METADATAHELPER_H_
