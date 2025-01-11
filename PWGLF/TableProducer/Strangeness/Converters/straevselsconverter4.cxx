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
#include "ITStracking/Vertexer.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;

// Converts Stra Event selections from 000 to 001
struct straevselsconverter4 {
  Produces<aod::StraEvSels_004> straEvSels_004;

  void process(aod::StraEvSels_003 const& straEvSels_003)
  {
    for (auto& values : straEvSels_003) {
      straEvSels_004(values.sel8(),
                     values.selection_raw(),
                     values.multFT0A(),
                     values.multFT0C(),
                     values.multFT0A(),
                     0 /*dummy FDDA value*/,
                     0 /*dummy FDDC value*/,
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
                     values.trackOccupancyInTimeRange(),
                     0 /*dummy occupancy value*/,
                     values.gapSide(),
                     values.totalFT0AmplitudeA(),
                     values.totalFT0AmplitudeC(),
                     values.totalFV0AmplitudeA(),
                     values.totalFDDAmplitudeA(),
                     values.totalFDDAmplitudeC(),
                     values.energyCommonZNA(),
                     values.energyCommonZNC(),
                     values.flags(),
                     0 /*dummy Alias value*/);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselsconverter4>(cfgc)};
}
