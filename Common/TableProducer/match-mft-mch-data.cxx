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

/// \file match-mft-mch-data.cxx
/// \brief Data entry point for the shared MFT MCH workflow.
/// \author Satoshi Yano <satoshi.yano@cern.ch>

// o2-linter: disable=name/workflow-file (Thin entry point delegates to the shared workflow factory.)
// o2-linter: disable=name/file-cpp (Preserve the established matching task source path.)

#include "Common/Core/MftMchMatching.h"

void customize(std::vector<o2::framework::ConfigParamSpec>& options)
{
  options.push_back({"matching-output-mode", o2::framework::VariantType::String, "auto", {"auto reads matching-output.cfgOutputMode from local JSON; tree, sparse or both overrides it"}});
}

#include <Framework/runDataProcessing.h>

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  return o2::analysis::mftmch::makeMatchingWorkflow(context);
}
