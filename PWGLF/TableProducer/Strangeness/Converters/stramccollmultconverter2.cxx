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

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct stramccollmultconverter2 {
  Produces<aod::StraMCCollMults_002> straMCCollMults_002;

  void process(aod::StraMCCollMults_001 const& straMCcolls)
  {
    straMCCollMults_002.reserve(straMCcolls.size());
    for (auto& straMCcoll : straMCcolls) {
      straMCCollMults_002(straMCcoll.multMCFT0A(),
                          straMCcoll.multMCFT0C(),
                          -1., // dummy value multMCFV0A
                          -1., // dummy value multMCFDDA
                          -1., // dummy value multMCFDDC
                          straMCcoll.multMCNParticlesEta05(),
                          straMCcoll.multMCNParticlesEta08(),
                          straMCcoll.multMCNParticlesEta10(),
                          straMCcoll.totalMultMCParticles());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stramccollmultconverter2>(cfgc)};
}
