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
/// \author Junlee Kim (jikim1290@gmail.com)
/// \since November 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct lumiTask {
  HistogramRegistry histos{"histos", {
                                       {"vertexx", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},                                //
                                       {"vertexy", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},                                //
                                       {"timestamp", "", {HistType::kTH1F, {{1000, 0, 5e7, "t"}}}},                             //
                                       {"vertexx_timestamp", "", {HistType::kTH2F, {{1000, 0, 5e7, "t"}, {1000, -1, 1, "x"}}}}, //
                                       {"vertexy_timestamp", "", {HistType::kTH2F, {{1000, 0, 5e7, "t"}, {1000, -1, 1, "y"}}}}  //
                                     }};
  uint64_t first_time = 1530314294062; // to be updated

  void process(aod::Collision const& collision, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    histos.fill(HIST("vertexx"), collision.posX());
    histos.fill(HIST("vertexy"), collision.posY());
    histos.fill(HIST("timestamp"), bc.timestamp() - first_time);
    //		LOGP(info, "Got timestamp {}", bc.timestamp());
    histos.fill(HIST("vertexx_timestamp"), bc.timestamp() - first_time, collision.posX());
    histos.fill(HIST("vertexy_timestamp"), bc.timestamp() - first_time, collision.posY());
  } // need selections
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lumiTask>(cfgc)};
}
