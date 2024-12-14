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
//
// This code calculates output histograms for centrality calibration
// as well as vertex-Z dependencies of raw variables (either for calibration
// of vtx-Z dependencies or for the calibration of those).
//
// This task is not strictly necessary in a typical analysis workflow,
// except for centrality calibration! The necessary task is the multiplicity
// tables.
//
// Comments, suggestions, questions? Please write to:
// - victor.gonzalez@cern.ch
// - david.dobrigkeit.chinellato@cern.ch
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

struct OTFV0Qa {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};

    // Base histograms
    histos.add("hEventCounter", "Event counter", kTH1D, {axisEvent});
  }

  void process(aod::Collision const& col)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OTFV0Qa>(cfgc)};
}
