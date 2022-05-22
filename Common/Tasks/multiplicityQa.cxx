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
  Configurable<bool> INELgtZERO{"INELgtZERO", 1, "0 - no, 1 - yes"};

  Configurable<int> NBinsMultFV0{"NBinsMultFV0", 1000, "N bins FV0"};
  Configurable<int> NBinsMultFT0{"NBinsMultFT0", 1000, "N bins FT0"};
  Configurable<int> NBinsMultFDD{"NBinsMultFDD", 1000, "N bins FDD"};
  Configurable<int> NBinsMultNTracks{"NBinsMultNTracks", 1000, "N bins Ntracks"};
  Configurable<float> MaxMultFV0{"MaxMultFV0", 20000, "Max FV0 signal"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 10000, "Max FT0 signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 10000, "Max FDD signal"};
  Configurable<float> MaxMultNTracks{"MaxMultNTracks", 1000, "Max Ntracks"};
  Configurable<int> NBinsVertexZ{"NBinsVertexZ", 400, "max vertex Z (cm)"};

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};
    const AxisSpec axisMultFV0{(int)NBinsMultFV0, 0, MaxMultFV0, "FV0 total amplitude"};
    const AxisSpec axisMultFT0{(int)NBinsMultFT0, 0, MaxMultFT0, "FT0 total amplitude"};
    const AxisSpec axisMultFDD{(int)NBinsMultFDD, 0, MaxMultFDD, "FDD total amplitude"};
    const AxisSpec axisMultNTracks{(int)NBinsMultNTracks, 0, MaxMultNTracks, "Track counter"};
    const AxisSpec axisVertexZ{(int)NBinsVertexZ, -20, 20, "Vertex Z (cm)"};

    //Base histograms
    histos.add("multiplicityQa/hEventCounter", "Event counter", kTH1D, {axisEvent});
    histos.add("multiplicityQa/hRawFV0", "Raw FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hRawFT0", "Raw FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hRawFDD", "Raw FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hRawNTracksPV", "Raw NTracks", kTH1D, {axisMultNTracks});
    histos.add("multiplicityQa/hZeqFV0", "vtx-z eq FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hZeqFT0", "vtx-z eq FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hZeqFDD", "vtx-z eq FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hZeqNTracksPV", "vtx-z eq NTracks", kTH1D, {axisMultNTracks});

    //Vertex-Z profiles for vertex-Z dependency estimate
    histos.add("multiplicityQa/hVtxZFV0A", "Av FV0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0A", "Av FT0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0C", "Av FT0C vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDA", "Av FDDA vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDC", "Av FDDC vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZNTracksPV", "Av NTracks vs vertex Z", kTProfile, {axisVertexZ});
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>::iterator const& col)
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
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    histos.fill(HIST("multiplicityQa/hEventCounter"), 2.5);

    //Vertex-Z dependencies, necessary for CCDB objects
    histos.fill(HIST("multiplicityQa/hVtxZFV0A"), col.posZ(), col.multFV0A());
    histos.fill(HIST("multiplicityQa/hVtxZFT0A"), col.posZ(), col.multFT0A());
    histos.fill(HIST("multiplicityQa/hVtxZFT0C"), col.posZ(), col.multFT0C());
    histos.fill(HIST("multiplicityQa/hVtxZFDDA"), col.posZ(), col.multFDDA());
    histos.fill(HIST("multiplicityQa/hVtxZFDDC"), col.posZ(), col.multFDDC());
    histos.fill(HIST("multiplicityQa/hVtxZNTracksPV"), col.posZ(), col.multNTracksPV());

    if (fabs(col.posZ()) > vtxZsel) {
      return;
    }

    histos.fill(HIST("multiplicityQa/hEventCounter"), 3.5);

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFV0M=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFT0M=%5.0f multFDDA=%5.0f multFDDC=%5.0f", col.multFV0A(), col.multFV0C(), col.multFV0M(), col.multFT0A(), col.multFT0C(), col.multFT0M(), col.multFDDA(), col.multFDDC());

    //Raw multiplicities
    histos.fill(HIST("multiplicityQa/hRawFV0"), col.multFV0A());
    histos.fill(HIST("multiplicityQa/hRawFT0"), col.multFT0M());
    histos.fill(HIST("multiplicityQa/hRawFDD"), col.multFDDM());
    histos.fill(HIST("multiplicityQa/hRawNTracksPV"), col.multNTracksPV());

    //vertex-Z corrected - FIXME
    histos.fill(HIST("multiplicityQa/hZeqFV0"), col.multZeqFV0A());
    histos.fill(HIST("multiplicityQa/hZeqFT0"), col.multZeqFT0A() + col.multZeqFT0C());
    histos.fill(HIST("multiplicityQa/hZeqFDD"), col.multZeqFDDA() + col.multZeqFDDC());
    histos.fill(HIST("multiplicityQa/hZeqNTracksPV"), col.multZeqNTracksPV());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityQa>(cfgc)};
}
