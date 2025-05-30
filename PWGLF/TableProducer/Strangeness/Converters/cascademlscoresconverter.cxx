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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct cascademlscoresconverter {
  Produces<aod::CascXiMLScores> xiMLSelections;
  Produces<aod::CascOmMLScores> omegaMLSelections;

  void process(aod::CascCores const& casccores)
  {
    for (int64_t i = 0; i < casccores.size(); ++i) {
      xiMLSelections(-1);
      omegaMLSelections(-1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascademlscoresconverter>(cfgc)};
}
