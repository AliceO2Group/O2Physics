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
/// \file   TableHelper.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base of utilities to build advanced tasks
///

#include "Common/Core/TableHelper.h"

#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>

#include <string>

/// Function to print the table required in the full workflow
/// @param initContext initContext of the init function
void o2::common::core::printTablesInWorkflow(o2::framework::InitContext& initContext)
{
  const auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  for (auto const& device : workflows.devices) {
    for (auto const& input : device.inputs) {
      LOG(info) << "Table: " << input.matcher.binding << " in device: " << device.name;
    }
  }
}

/// Function to check if a table is required in a workflow
/// @param initContext initContext of the init function
/// @param table name of the table to check for
bool o2::common::core::isTableRequiredInWorkflow(o2::framework::InitContext& initContext, const std::string& table)
{
  LOG(debug) << "Checking if table " << table << " is needed";
  bool tableNeeded = false;
  const auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  for (auto const& device : workflows.devices) {
    for (auto const& input : device.inputs) {
      if (input.matcher.binding == table) {
        LOG(debug) << "Table: " << input.matcher.binding << " is needed in device: " << device.name;
        tableNeeded = true;
      }
    }
  }
  return tableNeeded;
}

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag bool value of flag to set, if the given value is true it will be kept, disregarding the table usage in the workflow.
void o2::common::core::enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, bool& flag)
{
  if (flag) {
    LOG(info) << "Table enabled: " + table;
    return;
  }
  if (isTableRequiredInWorkflow(initContext, table)) {
    flag = true;
    LOG(info) << "Auto-enabling table: " + table;
    return;
  }
  LOG(info) << "Table disabled and not required: " + table;
}

/// Function to enable or disable a configurable flag, depending on the fact that a table is needed or not
/// @param initContext initContext of the init function
/// @param table name of the table to check for
/// @param flag int value of flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
void o2::common::core::enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, int& flag)
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
      return;
    }
    LOG(warn) << "Table disabled but required: " + table;
  } else {
    if (flag < 0) {
      flag = 0;
      LOG(info) << "Auto-disabling table: " + table;
      return;
    }
    LOG(info) << "Table disabled and not required: " + table;
  }
}
