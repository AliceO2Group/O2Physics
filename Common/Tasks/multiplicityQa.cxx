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
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

struct MultiplicityQa {
  //Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"sel", 7, "trigger: 7 - sel7, 8 - sel8"};
  Configurable<float> vtxZsel{"vtxZsel", 10, "max vertex Z (cm)"};
  
  Configurable<int> NBinsMultV0{"NBinsMultV0", 1000, "N bins V0"};
  Configurable<int> NBinsMultT0{"NBinsMultT0", 1000, "N bins T0"};
  Configurable<int> NBinsMultFDD{"NBinsMultFDD", 1000, "N bins FDD"};
  Configurable<int> NBinsMultNTracks{"NBinsMultNTracks", 1000, "N bins Ntracks"};
  Configurable<float> MaxMultV0{"MaxMultV0", 20000, "Max V0 signal"};
  Configurable<float> MaxMultT0{"MaxMultT0", 10000, "Max T0 signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 10000, "Max FDD signal"};
  Configurable<float> MaxMultNTracks{"MaxMultNTracks", 1000, "Max Ntracks"};
  Configurable<int> NBinsVertexZ{"NBinsVertexZ", 400, "max vertex Z (cm)"};
  
  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0,10, "Event counter"};
    const AxisSpec axisMultV0{(int)NBinsMultV0, 0, MaxMultV0, "V0 total amplitude"};
    const AxisSpec axisMultT0{(int)NBinsMultT0, 0, MaxMultT0, "T0 total amplitude"};
    const AxisSpec axisMultFDD{(int)NBinsMultFDD, 0, MaxMultFDD, "FDD total amplitude"};
    const AxisSpec axisMultNTracks{(int)NBinsMultNTracks, 0, MaxMultNTracks, "Track counter"};
    const AxisSpec axisVertexZ{(int)NBinsVertexZ, -20,20, "Vertex Z (cm)"};
    
    //Base histograms
    histos.add("multiplicityQa/hEventCounter", "Event counter", kTH1D, {axisEvent});
    histos.add("multiplicityQa/hRawV0", "Raw V0", kTH1D, {axisMultV0});
    histos.add("multiplicityQa/hRawT0", "Raw T0", kTH1D, {axisMultT0});
    histos.add("multiplicityQa/hRawFDD", "Raw FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hRawNTracks", "Raw NTracks", kTH1D, {axisMultNTracks});
    histos.add("multiplicityQa/hZeqV0", "vtx-z eq V0", kTH1D, {axisMultV0});
    histos.add("multiplicityQa/hZeqT0", "vtx-z eq T0", kTH1D, {axisMultT0});
    histos.add("multiplicityQa/hZeqFDD", "vtx-z eq FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hZeqNTracks", "vtx-z eq NTracks", kTH1D, {axisMultNTracks});
    
    //Vertex-Z profiles for vertex-Z dependency estimate
    histos.add("multiplicityQa/hVtxZV0A", "Av V0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZT0A", "Av T0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZT0C", "Av T0C vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDA", "Av FDDA vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDC", "Av FDDC vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZNTracks", "Av NTracks vs vertex Z", kTProfile, {axisVertexZ});
  }
  
  
  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& col)
  {
    histos.fill(HIST("multiplicityQa/hEventCounter"), 0.5);
    if (selection == 7 && !col.sel7()) {
      return;
    }

    if (selection == 8 && !col.sel8()) {
      return;
    }
    if (selection != 7 && selection != 8) {
      LOGF(fatal, "Unknown selection type! Use `--sel 7` or `--sel 8`");
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 1.5);
    
    //Vertex-Z dependencies, necessary for CCDB objects
    histos.fill(HIST("multiplicityQa/hVtxZV0A"), col.posZ(), col.multV0A());
    histos.fill(HIST("multiplicityQa/hVtxZT0A"), col.posZ(), col.multT0A());
    histos.fill(HIST("multiplicityQa/hVtxZT0C"), col.posZ(), col.multT0C());
    histos.fill(HIST("multiplicityQa/hVtxZFDDA"), col.posZ(), col.multFDDA());
    histos.fill(HIST("multiplicityQa/hVtxZFDDC"), col.posZ(), col.multFDDC());
    histos.fill(HIST("multiplicityQa/hVtxZNTracks"), col.posZ(), col.multNContribs());
    
    if (fabs(col.posZ()) > vtxZsel){
      return;
    }
    
    histos.fill(HIST("multiplicityQa/hEventCounter"), 2.5);

    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multV0M=%5.0f multT0A=%5.0f multT0C=%5.0f multT0M=%5.0f multFDDA=%5.0f multFDDC=%5.0f", col.multV0A(), col.multV0C(), col.multV0M(), col.multT0A(), col.multT0C(), col.multT0M(), col.multFDDA(), col.multFDDC());

    //Raw multiplicities
    histos.fill(HIST("multiplicityQa/hRawV0"), col.multV0A());
    histos.fill(HIST("multiplicityQa/hRawT0"), col.multT0M());
    histos.fill(HIST("multiplicityQa/hRawFDD"), col.multFDD());
    histos.fill(HIST("multiplicityQa/hRawNTracks"), col.multNContribs());
    
    //vertex-Z corrected - FIXME
    histos.fill(HIST("multiplicityQa/hZeqV0"), col.multV0A());
    histos.fill(HIST("multiplicityQa/hZeqT0"), col.multZeqT0A()+col.multZeqT0C());
    histos.fill(HIST("multiplicityQa/hZeqFDD"), col.multZeqFDDA()+col.multZeqFDDC());
    histos.fill(HIST("multiplicityQa/hZeqNTracks"), col.multZeqNContribs());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityQa>(cfgc)};
}
