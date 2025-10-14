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

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct Run2TracksExtraConverter {
  Produces<aod::Run2TrackExtras_001> Run2TrackExtras_001;
  void process(aod::Run2TrackExtras_000 const& Run2TrackExtras_000)
  {

    for (const auto& track0 : Run2TrackExtras_000) {
      Run2TrackExtras_001(track0.itsSignal(), 0);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Run2TracksExtraConverter>(cfgc)};
}
