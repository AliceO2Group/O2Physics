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
struct straevselsconverter2rawcents2 {
  Produces<aod::StraRawCents_004> straRawCents_004;

  void process(aod::StraEvSels_002 const& straEvSels_002)
  {
    for (auto& values : straEvSels_002) {
      straRawCents_004(values.multFT0A(),
                       values.multFT0C(),
                       values.multFT0A(),
                       values.multNTracksPVeta1(),
                       values.multPVTotalContributors(),
                       values.multNTracksGlobal(),
                       values.multNTracksITSTPC(),
                       values.multAllTracksTPCOnly(),
                       values.multAllTracksITSTPC(),
                       values.multZNA(),
                       values.multZNC(),
                       values.multZEM1(),
                       values.multZEM2(),
                       values.multZPA(),
                       values.multZPC(),
                       values.trackOccupancyInTimeRange());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselsconverter2rawcents2>(cfgc)};
}
