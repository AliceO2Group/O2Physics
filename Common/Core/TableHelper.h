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

#include <Framework/Array2D.h>
#include <Framework/ConfigParamSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>

#include <string>
#include <vector>

namespace o2::common::core
{

/// Function to print the table required in the full workflow
/// @param initContext initContext of the init function
void printTablesInWorkflow(o2::framework::InitContext& initContext);

/// Function to check if a table is required in a workflow
/// @param initContext initContext of the init function
/// @param table name of the table to check for
bool isTableRequiredInWorkflow(o2::framework::InitContext& initContext, const std::string& table);

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag bool value of flag to set, if the given value is true it will be kept, disregarding the table usage in the workflow.
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, bool& flag);

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag int value of flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, int& flag);

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag configurable flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
template <typename FlagType>
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, FlagType& flag)
{
  enableFlagIfTableRequired(initContext, table, flag.value);
}

/// Function to check for a specific configurable from another task in the current workflow and fetch its value. Useful for tasks that need to know the value of a configurable in another task.
/// @param initContext initContext of the init function
/// @param taskName name of the task to check for
/// @param optName name of the option to check for
/// @param value value of the option to set
/// @param verbose if true, print debug messages
template <typename ValueType>
bool getTaskOptionValue(o2::framework::InitContext& initContext, const std::string& taskName, const std::string& optName, ValueType& value, const bool verbose = true)
{
  if (verbose) {
    LOG(info) << "Checking for option '" << optName << "' in task '" << taskName << "'";
  }
  const auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  int deviceCounter = 0;
  bool found = false;
  for (auto const& device : workflows.devices) {
    if (verbose) {
      LOG(info) << " Device " << deviceCounter++ << " " << device.name;
    }
    if (device.name == taskName) { // Found the mother task
      int optionCounter = 0;
      for (const o2::framework::ConfigParamSpec& option : device.options) {
        if (verbose) {
          LOG(info) << "  Option " << optionCounter++ << " " << option.name << " of type " << static_cast<int>(option.type) << " = '" << option.defaultValue.asString() << "'";
        }
        if (option.name == optName) {
          value = option.defaultValue.get<ValueType>();
          if (verbose) {
            if constexpr (!std::is_same_v<ValueType, o2::framework::LabeledArray<float>>) {
              LOG(info) << "   Found option '" << optName << " of type LabeledArray<float>";
            } else if constexpr (!std::is_same_v<ValueType, std::vector<std::string>>) {
              LOG(info) << "   Found option '" << optName << " of type vector<string>";
            } else {
              LOG(info) << "   Found option '" << optName << "' with value '" << value << "'";
            }
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

/// Function to check for a specific configurable from another task in the current workflow and fetch its value. Useful for tasks that need to know the value of a configurable in another task.
/// @param initContext initContext of the init function
/// @param taskName name of the task to check for
/// @param value Task configurable to inherit from (name and values are used)
/// @param verbose if true, print debug messages
template <typename ValueType>
bool getTaskOptionValue(o2::framework::InitContext& initContext, const std::string& taskName, ValueType& configurable, const bool verbose = true)
{
  return getTaskOptionValue(initContext, taskName, configurable.name, configurable.value, verbose);
}

} // namespace o2::common::core

using o2::common::core::enableFlagIfTableRequired;
using o2::common::core::getTaskOptionValue;
using o2::common::core::isTableRequiredInWorkflow;
using o2::common::core::printTablesInWorkflow;

#endif // COMMON_CORE_TABLEHELPER_H_
