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
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct stramccollisionconverter {
  Produces<aod::StraMCCollisions_001> straMCCollisions_001;

  void process(aod::StraMCCollisions_000 const& straMCcoll)
  {
    for (auto& mccollision : straMCcoll) {
      straMCCollisions_001(mccollision.posX(), mccollision.posY(), mccollision.posZ(),
                           mccollision.impactParameter(), 0.0f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stramccollisionconverter>(cfgc)};
}
