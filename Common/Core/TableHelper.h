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
/// \brief  Base of utilities to build advanced tasks
///

#ifndef COMMON_CORE_TABLEHELPER_H_
#define COMMON_CORE_TABLEHELPER_H_

#include <string>

#include "Framework/InitContext.h"
#include "Framework/RunningWorkflowInfo.h"

/// Function to check if a table is required in a workflow
/// @param initContext initContext of the init function
/// @param table name of the table to check for
bool isTableRequiredInWorkflow(o2::framework::InitContext& initContext, const std::string& table)
{
  auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  for (auto const& device : workflows.devices) {
    for (auto const& input : device.inputs) {
      if (input.matcher.binding == table) {
        return true;
      }
    }
  }
  return false;
}

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag int value of flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, int& flag)
{
  if (flag > 0) {
    flag = 1;
    LOG(info) << "Table enabled: " + table;
    return;
  }
  if (isTableRequiredInWorkflow(initContext, table)) {
    if (flag < 0) {
      flag = 1;
      LOG(info) << "Auto-enabling table: " + table;
    } else {
      LOG(info) << "Table disabled: " + table;
    }
  }
}

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag configurable flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
template <typename FlagType>
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, FlagType& flag)
{
  enableFlagIfTableRequired(initContext, table, flag.value);
}

/// Function to check for a specific configurable in the current workflow
template <typename ValueType>
bool getTaskOptionValue(o2::framework::InitContext& initContext, const std::string& taskName, const std::string& optName, ValueType& value, const bool verbose = true)
{
  auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  int deviceCounter = 0;
  bool found = false;
  for (auto const& device : workflows.devices) {
    if (verbose) {
      LOG(info) << " Device " << deviceCounter++ << " " << device.name;
    }
    if (device.name == taskName) { // Found the mother task
      int optionCounter = 0;
      for (auto const& option : device.options) {
        if (verbose) {
          LOG(info) << "  Option " << optionCounter++ << " " << option.name << " = '" << option.defaultValue.asString() << "'";
        }
        if (option.name == optName) {
          value = option.defaultValue.get<ValueType>();
          if (verbose) {
            LOG(info) << "   Found option '" << optName << "' with value '" << value << "'";
            found = true;
          } else {
            return true;
          }
        }
      }
      return found;
    }
  }
  return false;
}

#endif // COMMON_CORE_TABLEHELPER_H_
