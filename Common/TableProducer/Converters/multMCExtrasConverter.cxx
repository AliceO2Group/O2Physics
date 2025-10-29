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

#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct MultMCExtrasConverter {
  Produces<aod::MultMCExtras_001> multMCExtras_001;
  void process(aod::MultMCExtras_000 const& multMCExtras_000)
  {
    for (const auto& r : multMCExtras_000) {
      multMCExtras_001(r.multMCFT0A(), r.multMCFT0C(), 0, 0, 0,
                       r.multMCNParticlesEta05(),
                       r.multMCNParticlesEta08(),
                       r.multMCNParticlesEta10(),
                       r.multMCPVz());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultMCExtrasConverter>(cfgc)};
}
