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
// Questions / comments / suggestions? Please write to:
//  - victor.gonzalez@cern.ch
//  - david.dobrigkeit.chinellato@cern.ch
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

struct MultiplicityQa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"sel", 7, "trigger: 0 - none, 7 - sel7, 8 - sel8"};
  Configurable<float> vertexZ{"vertexZ", 10, "max vertex Z position (cm)"};
  
  //Configurables related to binning of the final histograms
  Configurable<float> MaxMultV0{"MaxMultV0", 10000.f, "Maximum multiplicity in V0"};
  Configurable<float> MaxMultT0{"MaxMultT0", 10000.f, "Maximum multiplicity in T0"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 10000.f, "Maximum multiplicity in FDD"};
  Configurable<int> NBinsMultV0{"NBinsMultV0", 1000, "N bins V0"};
  Configurable<int> NBinsMultT0{"NBinsMultT0", 1000, "N bins T0"};
  Configurable<int> NBinsMultFDD{"NBinsMultFDD", 1000, "N bins FDD"};
  Configurable<int> NBinsVertexZ{"NBinsVertexZ", 150, "N bins vertex Z profiles"};
  

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};
    const AxisSpec axisMultV0{(int)NBinsMultV0, 0, MaxMultV0, "V0 total amplitude"};
    const AxisSpec axisMultT0{(int)NBinsMultT0, 0, MaxMultT0, "T0 total amplitude"};
    const AxisSpec axisMultFDD{(int)NBinsMultFDD, 0, MaxMultFDD, "FDD total amplitude"};
    const AxisSpec axisMultZ{(int)NBinsVertexZ, -15, 15, "vertex Z (cm)"};
    
    //Initialize basic raw histograms
    histos.add("multiplicityQa/hEventCounter", "Event Counter", kTH1D, {axisEvent});
    histos.add("multiplicityQa/hRawV0", "V0 total amplitude (raw)", kTH1D, {axisMultV0});
    histos.add("multiplicityQa/hRawT0", "T0 total amplitude (raw)", kTH1D, {axisMultT0});
    histos.add("multiplicityQa/hRawFDD", "FDD total amplitude (raw)", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hZEqV0", "V0 total amplitude (vtx-z-eq)", kTH1D, {axisMultV0});
    histos.add("multiplicityQa/hZEqT0", "T0 total amplitude (vtx-z-eq)", kTH1D, {axisMultT0});
    histos.add("multiplicityQa/hZEqFDD", "FDD total amplitude (vtx-z-eq)", kTH1D, {axisMultFDD});
    
    //Initialize basic vertex-Z histograms
    histos.add("multiplicityQa/hVertexZV0A", "V0A vertex Z profile", kTProfile, {axisMultZ});
    histos.add("multiplicityQa/hVertexZT0A", "T0A vertex Z profile", kTProfile, {axisMultZ});
    histos.add("multiplicityQa/hVertexZT0C", "T0C vertex Z profile", kTProfile, {axisMultZ});
    histos.add("multiplicityQa/hVertexZFDDA", "FDDA vertex Z profile", kTProfile, {axisMultZ});
    histos.add("multiplicityQa/hVertexZFDDC", "FDDC vertex Z profile", kTProfile, {axisMultZ});
  }
  
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& col)
  {
    histos.fill(HIST("multiplicityQa/hEventCounter"), 0.5);
    if (!isMC && !col.alias()[kINT7]) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 1.5);
    if (selection == 7 && !col.sel7()) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 2.5);
    if (selection == 8 && !col.sel8()) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 3.5);
    if (selection != 7 && selection != 8 && selection != 0) {
      LOGF(fatal, "Unknown selection type! Use `--sel 7` or `--sel 8` or `--sel 0` (no selection)");
    }
    
    histos.fill(HIST("multiplicityQa/hVertexZV0A"), col.posZ(), col.multFDD());
    histos.fill(HIST("multiplicityQa/hVertexZT0A"), col.posZ(), col.multT0A());
    histos.fill(HIST("multiplicityQa/hVertexZT0C"), col.posZ(), col.multT0C());
    histos.fill(HIST("multiplicityQa/hVertexZFDDA"), col.posZ(), col.multFDDA());
    histos.fill(HIST("multiplicityQa/hVertexZFDDC"), col.posZ(), col.multFDDC());
    
    if (TMath::Abs(col.posZ())>vertexZ) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 4.5);
    
    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multV0M=%5.0f multT0A=%5.0f multT0C=%5.0f multT0M=%5.0f", col.multV0A(), col.multV0C(), col.multV0M(), col.multT0A(), col.multT0C(), col.multT0M());

    //Raw multiplicities
    histos.fill(HIST("multiplicityQa/hRawV0"), col.multV0M());
    histos.fill(HIST("multiplicityQa/hRawT0"), col.multT0M());
    histos.fill(HIST("multiplicityQa/hRawFDD"), col.multFDD());
    
    //Vertex-Z equalized multiplicities - FIXME
    histos.fill(HIST("multiplicityQa/hZEqV0"), col.multV0M());
    histos.fill(HIST("multiplicityQa/hZEqT0"), col.multT0M());
    histos.fill(HIST("multiplicityQa/hZEqFDD"), col.multFDD());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityQa>(cfgc)};
}
