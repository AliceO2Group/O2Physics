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
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

//STEP 0
//This is an empty analysis skeleton: the starting point! 
struct taskskeleton {
  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> etaHistogram{TH1F("etaHistogram", "etaHistogram", 200, -1., +1)};

  void process(aod::TracksIU const& tracks)
  {
    for (auto& track : tracks) {
      etaHistogram->Fill(track.eta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskskeleton>(cfgc)
  };
}
