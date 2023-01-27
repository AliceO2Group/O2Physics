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

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct MultiplicityQa {
  //Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"sel", 7, "trigger: 7 - sel7, 8 - sel8"};
  Configurable<float> vtxZsel{"vtxZsel", 10, "max vertex Z (cm)"};
  Configurable<bool> INELgtZERO{"INELgtZERO", 1, "0 - no, 1 - yes"};

  Configurable<int> NBinsMultFV0{"NBinsMultFV0", 1000, "N bins FV0"};
  Configurable<int> NBinsMultFT0{"NBinsMultFT0", 1000, "N bins FT0"};
  Configurable<int> NBinsMultFT0A{"NBinsMultFT0A", 1000, "N bins FT0A"};
  Configurable<int> NBinsMultFT0C{"NBinsMultFT0C", 1000, "N bins FT0C"};
  Configurable<int> NBinsMultFDD{"NBinsMultFDD", 1000, "N bins FDD"};
  Configurable<int> NBinsMultNTracks{"NBinsMultNTracks", 1000, "N bins Ntracks"};

  Configurable<int> NBinsMultFV02d{"NBinsMultFV02d", 100, "N bins FV0 in 2D"};
  Configurable<int> NBinsMultFT02d{"NBinsMultFT02d", 100, "N bins FT0 in 2D"};
  Configurable<int> NBinsMultFT0A2d{"NBinsMultFT0A2d", 100, "N bins FT0A in 2D"};
  Configurable<int> NBinsMultFT0C2d{"NBinsMultFT0C2d", 100, "N bins FT0C in 2D"};
  Configurable<int> NBinsMultFDD2d{"NBinsMultFDD2d", 100, "N bins FDD in 2D"};
  Configurable<int> NBinsMultNTracks2d{"NBinsMultNTracks2d", 100, "N bins Ntracks in 2D"};
  Configurable<int> NBinsNContributors{"NBinsNContributors", 100, "N bins Ntracks in 2D"};

  Configurable<float> MaxMultFV0{"MaxMultFV0", 20000, "Max FV0 signal"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 10000, "Max FT0 signal"};
  Configurable<float> MaxMultFT0A{"MaxMultFT0A", 10000, "Max FT0A signal"};
  Configurable<float> MaxMultFT0C{"MaxMultFT0C", 10000, "Max FT0C signal"};
  Configurable<float> MaxMultFDD{"MaxMultFDD", 10000, "Max FDD signal"};
  Configurable<float> MaxMultNTracks{"MaxMultNTracks", 1000, "Max Ntracks"};
  Configurable<float> MaxNContributors{"MaxNContributors", 200, "Max NContributors"};
  Configurable<int> NBinsVertexZ{"NBinsVertexZ", 400, "max vertex Z (cm)"};
  Configurable<bool> useZeqInProfiles{"useZeqInProfiles", true, "use Z-equalized signals in midrap Nch profiles"};

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};
    const AxisSpec axisMultFV0{(int)NBinsMultFV0, 0, MaxMultFV0, "FV0 total amplitude"};
    const AxisSpec axisMultFT0{(int)NBinsMultFT0, 0, MaxMultFT0, "FT0 total amplitude"};
    const AxisSpec axisMultFT0A{(int)NBinsMultFT0A, 0, MaxMultFT0A, "FT0A total amplitude"};
    const AxisSpec axisMultFT0C{(int)NBinsMultFT0C, 0, MaxMultFT0C, "FT0C total amplitude"};
    const AxisSpec axisMultFDD{(int)NBinsMultFDD, 0, MaxMultFDD, "FDD total amplitude"};
    const AxisSpec axisMultNTracks{(int)NBinsMultNTracks, 0, MaxMultNTracks, "Track counter"};

    const AxisSpec axisMultFV02d{(int)NBinsMultFV02d, 0, MaxMultFV0, "FV0 total amplitude"};
    const AxisSpec axisMultFT02d{(int)NBinsMultFT02d, 0, MaxMultFT0, "FT0 total amplitude"};
    const AxisSpec axisMultFT0A2d{(int)NBinsMultFT0A2d, 0, MaxMultFT0A, "FT0A total amplitude"};
    const AxisSpec axisMultFT0C2d{(int)NBinsMultFT0C2d, 0, MaxMultFT0C, "FT0C total amplitude"};
    const AxisSpec axisMultFDD2d{(int)NBinsMultFDD2d, 0, MaxMultFDD, "FDD total amplitude"};
    const AxisSpec axisMultNTracks2d{(int)NBinsMultNTracks2d, 0, MaxMultNTracks, "Track counter"};

    const AxisSpec axisVertexZ{(int)NBinsVertexZ, -20, 20, "Vertex Z (cm)"};
    const AxisSpec axisContributorsTRD{(int)NBinsNContributors, -0.5f, MaxNContributors - 0.5f, "N_{contribs}^{TRD}"};
    const AxisSpec axisContributorsTOF{(int)NBinsNContributors, -0.5f, MaxNContributors - 0.5f, "N_{contribs}^{TOF}"};

    //Base histograms
    histos.add("multiplicityQa/hEventCounter", "Event counter", kTH1D, {axisEvent});
    histos.add("multiplicityQa/hRawFV0", "Raw FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hRawFT0", "Raw FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hRawFT0A", "Raw FT0A", kTH1D, {axisMultFT0A});
    histos.add("multiplicityQa/hRawFT0C", "Raw FT0C", kTH1D, {axisMultFT0C});
    histos.add("multiplicityQa/hRawFDD", "Raw FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hRawNTracksPV", "Raw NTracks", kTH1D, {axisMultNTracks});
    histos.add("multiplicityQa/hZeqFV0", "vtx-z eq FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hZeqFT0", "vtx-z eq FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hZeqFT0A", "vtx-z eq FT0A", kTH1D, {axisMultFT0A});
    histos.add("multiplicityQa/hZeqFT0C", "vtx-z eq FT0C", kTH1D, {axisMultFT0C});
    histos.add("multiplicityQa/hZeqFDD", "vtx-z eq FDD", kTH1D, {axisMultFDD});
    histos.add("multiplicityQa/hZeqNTracksPV", "vtx-z eq NTracks", kTH1D, {axisMultNTracks});

    // Per BC
    histos.add("multiplicityQa/hPerBCRawFV0", "Raw FV0", kTH1D, {axisMultFV0});
    histos.add("multiplicityQa/hPerBCRawFT0", "Raw FT0", kTH1D, {axisMultFT0});
    histos.add("multiplicityQa/hPerBCRawFT0A", "Raw FT0A", kTH1D, {axisMultFT0A});
    histos.add("multiplicityQa/hPerBCRawFT0C", "Raw FT0C", kTH1D, {axisMultFT0C});
    histos.add("multiplicityQa/hPerBCRawFDD", "Raw FDD", kTH1D, {axisMultFDD});

    // Vertex-Z profiles for vertex-Z dependency estimate
    histos.add("multiplicityQa/hVtxZFV0A", "Av FV0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0A", "Av FT0A vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFT0C", "Av FT0C vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDA", "Av FDDA vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZFDDC", "Av FDDC vs vertex Z", kTProfile, {axisVertexZ});
    histos.add("multiplicityQa/hVtxZNTracksPV", "Av NTracks vs vertex Z", kTProfile, {axisVertexZ});

    // profiles of track contributors
    histos.add("multiplicityQa/hNchProfileFV0", "FV0", kTH2F, {axisMultFV02d, axisMultNTracks2d});
    histos.add("multiplicityQa/hNchProfileFT0", "FT0", kTH2F, {axisMultFT02d, axisMultNTracks2d});
    histos.add("multiplicityQa/hNchProfileFT0A", "FT0A", kTH2F, {axisMultFT0A2d, axisMultNTracks2d});
    histos.add("multiplicityQa/hNchProfileFT0C", "FT0C", kTH2F, {axisMultFT0C2d, axisMultNTracks2d});
    histos.add("multiplicityQa/hNchProfileFDD", "FDD", kTH2F, {axisMultFDD2d, axisMultNTracks2d});

    // Contributors correlation
    histos.add("h2dNContribCorrAll", "h2dNContribCorrAll", kTH2D, {axisContributorsTRD, axisContributorsTOF});
  }

  void processCollisions(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>::iterator const& col)
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
    histos.fill(HIST("multiplicityQa/hRawFT0A"), col.multFT0A());
    histos.fill(HIST("multiplicityQa/hRawFT0C"), col.multFT0C());
    histos.fill(HIST("multiplicityQa/hRawFDD"), col.multFDDM());
    histos.fill(HIST("multiplicityQa/hRawNTracksPV"), col.multNTracksPV());

    //vertex-Z corrected - FIXME
    histos.fill(HIST("multiplicityQa/hZeqFV0"), col.multZeqFV0A());
    histos.fill(HIST("multiplicityQa/hZeqFT0"), col.multZeqFT0A() + col.multZeqFT0C());
    histos.fill(HIST("multiplicityQa/hZeqFT0A"), col.multZeqFT0A());
    histos.fill(HIST("multiplicityQa/hZeqFT0C"), col.multZeqFT0C());
    histos.fill(HIST("multiplicityQa/hZeqFDD"), col.multZeqFDDA() + col.multZeqFDDC());
    histos.fill(HIST("multiplicityQa/hZeqNTracksPV"), col.multZeqNTracksPV());

    // Profiles
    if (useZeqInProfiles) {
      histos.fill(HIST("multiplicityQa/hNchProfileFV0"), col.multZeqFV0A(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0"), col.multZeqFT0A() + col.multZeqFT0C(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0A"), col.multZeqFT0A(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0C"), col.multZeqFT0C(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFDD"), col.multZeqFDDA() + col.multZeqFDDC(), col.multZeqNTracksPV());
    } else {
      histos.fill(HIST("multiplicityQa/hNchProfileFV0"), col.multFV0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0"), col.multFT0A() + col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0A"), col.multFT0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFT0C"), col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/hNchProfileFDD"), col.multFDDA() + col.multFDDC(), col.multNTracksPV());
    }
  }
  PROCESS_SWITCH(MultiplicityQa, processCollisions, "per-collision analysis", true);

  void processBCs(BCsWithRun3Matchings::iterator const& bc,
                  aod::FV0As const&,
                  aod::FT0s const&,
                  aod::FDDs const&)
  {
    float multFV0A = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;

    if (bc.has_ft0()) {
      auto ft0 = bc.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    }

    if (bc.has_fdd()) {
      auto fdd = bc.fdd();
      for (auto amplitude : fdd.chargeA()) {
        multFDDA += amplitude;
      }
      for (auto amplitude : fdd.chargeC()) {
        multFDDC += amplitude;
      }
    }
    // using FV0 row index from event selection task
    if (bc.has_fv0a()) {
      auto fv0a = bc.fv0a();
      for (auto amplitude : fv0a.amplitude()) {
        multFV0A += amplitude;
      }
    }

    // Raw multiplicities per BC
    histos.fill(HIST("multiplicityQa/hPerBCRawFV0"), multFV0A);
    histos.fill(HIST("multiplicityQa/hPerBCRawFT0"), multFT0A + multFT0C);
    histos.fill(HIST("multiplicityQa/hPerBCRawFT0A"), multFT0A);
    histos.fill(HIST("multiplicityQa/hPerBCRawFT0C"), multFT0C);
    histos.fill(HIST("multiplicityQa/hPerBCRawFDD"), multFDDA + multFDDC);
  }
  PROCESS_SWITCH(MultiplicityQa, processBCs, "per-BC analysis", true);

  void processCollisionsPVChecks(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>::iterator const& col, soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks)
  {
    if (selection == 7 && !col.sel7()) {
      return;
    }

    if (selection == 8 && !col.sel8()) {
      return;
    }
    if (selection != 7 && selection != 8) {
      LOGF(fatal, "Unknown selection type! Use `--sel 7` or `--sel 8`");
    }

    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (fabs(col.posZ()) > vtxZsel) {
      return;
    }

    long NcontribsTOF = 0;
    long NcontribsTRD = 0;
    for (auto& track : tracks) {
      if (track.isPVContributor()) {
        if (track.hasTRD())
          NcontribsTRD++;
        if (track.hasTOF())
          NcontribsTOF++;
      }
    }
    histos.fill(HIST("h2dNContribCorrAll"), NcontribsTRD, NcontribsTOF);
  }
  PROCESS_SWITCH(MultiplicityQa, processCollisionsPVChecks, "do PV contributors check", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityQa>(cfgc)};
}
