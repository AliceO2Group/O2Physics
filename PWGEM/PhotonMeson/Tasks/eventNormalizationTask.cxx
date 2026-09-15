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

/// \file eventNormalizationTask.cxx
/// \brief task to produce histogram with event counter information
/// \author Joshua Konig, joshua.konig@cern.ch

#include "PWGEM/PhotonMeson/DataModel/EventTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1D.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct EventNormalizationTask {

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext& /*context*/)
  {
    // Names derived from EventAcceptance bits in EventTables.h
    std::vector<std::string> vecLabelNames = {
      "all",
      "has MC coll",
      "good zVtx",
      "is FT0AND",
      "No TFB",
      "ITS ROFB",
      "Same bunch Pileup",
      "ZVtxFT0PV",
      "No coll in time range",
      "good track occi",
      "Good FT0 occupancy",
      "kTVXinEMC",
      "Good centrality",
      "Good RCT",
      "Good Sel8"};
    fRegistry.add("EventNormalization/hEventCounter", "Event Counter;category;#it{N}_{evt}", kTH1D, {{static_cast<int>(vecLabelNames.size()), -0.5, -0.5 + static_cast<int>(vecLabelNames.size())}}, false);
    for (size_t i = 0; i < vecLabelNames.size(); ++i) {
      fRegistry.get<TH1>(HIST("EventNormalization/hEventCounter"))->GetXaxis()->SetBinLabel(i + 1, vecLabelNames[i].c_str());
    }
  }

  // Process function only takes the entries and fills them into the histogram. This is done per DF and not per collision
  void processNorm(o2::aod::PMEvSelBit const& evSelBit)
  {

    auto vecEvSelBits = evSelBit.eventSelectionBit();
    for (size_t i = 0; i < vecEvSelBits.size(); ++i) {
      fRegistry.fill(HIST("EventNormalization/hEventCounter"), i, vecEvSelBits[i]);
    }
  }

  PROCESS_SWITCH(EventNormalizationTask, processNorm, "run event normalization task", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<EventNormalizationTask>(context, TaskName{"event-normalization-task"})};
}
