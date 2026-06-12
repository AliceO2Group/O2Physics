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

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

// Converts Stra Event selections from 005 to 006
struct straevselsconverter6 {
  Produces<aod::StraEvSels_006> straEvSels_006;

  void process(aod::StraEvSels_005 const& straEvSels_005)
  {
    for (auto& values : straEvSels_005) {
      straEvSels_006(values.sel8(),
                     values.selection_raw(),
                     values.multFT0A(),
                     values.multFT0C(),
                     values.multFT0A(),
                     values.multFDDA(),
                     values.multFDDC(),
                     values.multNTracksPVeta1(),
                     values.multPVTotalContributors(),
                     values.multNTracksGlobal(),
                     values.flags(),
                     values.alias_raw(),
                     values.rct_raw());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselsconverter6>(cfgc)};
}
