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

#ifndef O2PHYSICS_COMMON_CORE_TABLEHELPER_H_
#define O2PHYSICS_COMMON_CORE_TABLEHELPER_H_

#include "Framework/InitContext.h"
#include "Framework/RunningWorkflowInfo.h"
#include <string>

/// Function to check if a table is required in a workflow
/// @param initContext initContext of the init function
/// @param table name of the table to check for
bool isTableRequiredInWorkflow(o2::framework::InitContext& initContext, const std::string& table)
{
  auto& workflows = initContext.services().get<o2::framework::RunningWorkflowInfo const>();
  for (auto device : workflows.devices) {
    for (auto input : device.inputs) {
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
/// @param flag configurable flag to set, only if initially set to -1. Initial values of 0 or 1 will be kept disregarding the table usage in the workflow.
template <typename FlagType>
void enableFlagIfTableRequired(o2::framework::InitContext& initContext, const std::string& table, FlagType& flag)
{
  if (isTableRequiredInWorkflow(initContext, table)) {
    if (flag < 0) {
      flag.value = 1;
      LOG(info) << "Auto-enabling table: " + table;
    } else if (flag > 0) {
      flag.value = 1;
      LOG(info) << "Table enabled: " + table;
    } else {
      LOG(info) << "Table disabled: " + table;
    }
  }
}

#endif
