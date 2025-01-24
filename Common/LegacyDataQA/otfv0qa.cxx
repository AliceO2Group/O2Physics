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

  Configurable<float> maxGammaMassForXYplot{"maxGammaMassForXYplot", 0.1f, "Max photon mass for XY plot"};

  ConfigurableAxis axisNCandidates{"axisNCandidates", {500, 0, 500}, "Number of OTF v0s"};
  ConfigurableAxis axisPosition{"axisPosition", {1000, -100, 100}, "position (cm)"};
  ConfigurableAxis axisMass{"axisMass", {100, 0.0f, 1.0f}, "Mass (GeV/c2)"};

  void init(InitContext&)
  {
    const AxisSpec axisPVz{30, -15, 15, "Primary vertex Z (cm)"};

    // Base histograms
    histos.add("hPrimaryVertexZ", "Event counter", kTH1F, {axisPVz});
    histos.add("hCandidates", "Number of OTF V0s", kTH1F, {axisNCandidates});
    histos.add("hGammaMass", "mass distribution", kTH1F, {axisMass});
    histos.add("h2dPosition", "xy positions", kTH2F, {axisPosition, axisPosition});
  }

  void process(aod::Collision const& collision, aod::Run2OTFV0s const& v0s)
  {
    histos.fill(HIST("hPrimaryVertexZ"), collision.posZ());
    histos.fill(HIST("hCandidates"), v0s.size());
    for (auto const& v0 : v0s) {
      histos.fill(HIST("hGammaMass"), v0.mass());
      if (v0.mass() < maxGammaMassForXYplot) {
        histos.fill(HIST("h2dPosition"), v0.x(), v0.y());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OTFV0Qa>(cfgc)};
}
