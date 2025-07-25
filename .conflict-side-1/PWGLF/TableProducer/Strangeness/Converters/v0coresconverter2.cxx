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
struct v0coresconverter2 {
  Produces<aod::V0MCCores_002> v0MCCores_002;

  void process(aod::V0MCCores_001 const& v0MCCores_001)
  {
    for (auto& values : v0MCCores_001) {
      v0MCCores_002(0,
                    values.pdgCode(),
                    values.pdgCodeMother(),
                    values.pdgCodePositive(),
                    values.pdgCodeNegative(),
                    values.isPhysicalPrimary(),
                    values.xMC(),
                    values.yMC(),
                    values.zMC(),
                    values.pxPosMC(),
                    values.pyPosMC(),
                    values.pzPosMC(),
                    values.pxNegMC(),
                    values.pyNegMC(),
                    values.pzNegMC(),
                    values.pxPosMC() + values.pxNegMC(),
                    values.pyPosMC() + values.pyNegMC(),
                    values.pzPosMC() + values.pzNegMC());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0coresconverter2>(cfgc)};
}
