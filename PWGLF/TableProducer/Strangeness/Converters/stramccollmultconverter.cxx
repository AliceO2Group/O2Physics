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

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct stramccollmultconverter {
  Produces<aod::StraMCCollMults_001> straMCCollMults_001;

  void process(aod::StraMCCollMults_000 const& straMCcolls)
  {
    for (auto& straMCcoll : straMCcolls) {
      straMCCollMults_001(straMCcoll.multMCFT0A(),
                          straMCcoll.multMCFT0C(),
                          straMCcoll.multMCNParticlesEta05(),
                          straMCcoll.multMCNParticlesEta08(),
                          straMCcoll.multMCNParticlesEta10(),
                          -1 /* dummy value for totalMultMCParticles */);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stramccollmultconverter>(cfgc)};
}
