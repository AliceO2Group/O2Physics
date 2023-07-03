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
#include "PWGLF/DataModel/LFQATables.h"
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
  Configurable<bool> INELgtZERO{"INELgtZERO", true, "0 - no, 1 - yes"};
  Configurable<bool> do2Dplots{"do2Dplots", false, "0 - no, 1 - yes"};

  ConfigurableAxis axisMultFV0{"axisMultFV0", {10000, 0, 500000}, "FV0 amplitude"};
  ConfigurableAxis axisMultFT0{"axisMultFT0", {10000, 0, 40000}, "FT0 amplitude"};
  ConfigurableAxis axisMultFT0A{"axisMultFT0A", {10000, 0, 30000}, "FT0A amplitude"};
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {10000, 0, 10000}, "FT0C amplitude"};
  ConfigurableAxis axisMultFDD{"axisMultFDD", {1000, 0, 4000}, "FDD amplitude"};
  ConfigurableAxis axisMultNTracks{"axisMultNTracks", {500, 0, 500}, "N_{tracks}"};
  ConfigurableAxis axisVertexZ{"axisVertexZ", {60, -15, 15}, "Vertex z (cm)"};

  ConfigurableAxis axisContributors{"axisContributors", {100, -0.5f, 99.5f}, "Vertex z (cm)"};
  ConfigurableAxis axisNumberOfPVs{"axisNumberOfPVs", {10, -0.5f, 9.5f}, "Number of reconstructed PVs"};

  Configurable<bool> useZeqInProfiles{"useZeqInProfiles", true, "use Z-equalized signals in midrap Nch profiles"};

  void init(InitContext&)
  {
    const AxisSpec axisEvent{10, 0, 10, "Event counter"};

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

    // two-dimensional histograms
    if (do2Dplots) {
      histos.add("multiplicityQa/h2dNchVsFV0", "FV0", kTH2F, {axisMultFV0, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0", "FT0", kTH2F, {axisMultFT0, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0A", "FT0A", kTH2F, {axisMultFT0A, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFT0C", "FT0C", kTH2F, {axisMultFT0C, axisMultNTracks});
      histos.add("multiplicityQa/h2dNchVsFDD", "FDD", kTH2F, {axisMultFDD, axisMultNTracks});

      histos.add("multiplicityQa/h2dPVsVsFV0", "FV0", kTH2F, {axisMultFV0, axisNumberOfPVs});
      histos.add("multiplicityQa/h2dPVsVsFT0", "FT0", kTH2F, {axisMultFT0, axisNumberOfPVs});
      histos.add("multiplicityQa/h2dPVsVsFT0A", "FT0A", kTH2F, {axisMultFT0A, axisNumberOfPVs});
      histos.add("multiplicityQa/h2dPVsVsFT0C", "FT0C", kTH2F, {axisMultFT0C, axisNumberOfPVs});
      histos.add("multiplicityQa/h2dPVsVsFDD", "FDD", kTH2F, {axisMultFDD, axisNumberOfPVs});
    }

    // Contributors correlation
    histos.add("h2dNContribCorrAll", "h2dNContribCorrAll", kTH2D, {axisContributors, axisContributors});
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
    if (useZeqInProfiles && do2Dplots) {
      histos.fill(HIST("multiplicityQa/h2dNchVsFV0"), col.multZeqFV0A(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0"), col.multZeqFT0A() + col.multZeqFT0C(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0A"), col.multZeqFT0A(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0C"), col.multZeqFT0C(), col.multZeqNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFDD"), col.multZeqFDDA() + col.multZeqFDDC(), col.multZeqNTracksPV());
    } else {
      histos.fill(HIST("multiplicityQa/h2dNchVsFV0"), col.multFV0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0"), col.multFT0A() + col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0A"), col.multFT0A(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFT0C"), col.multFT0C(), col.multNTracksPV());
      histos.fill(HIST("multiplicityQa/h2dNchVsFDD"), col.multFDDA() + col.multFDDC(), col.multNTracksPV());
    }
  }

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

  void processMCCollisions(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs, aod::McCollisionLabels>::iterator const& col, soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions)
  {
    if (selection == 7 && !col.sel7()) {
      return;
    }

    if (selection == 8 && !col.sel8()) {
      return;
    }
    if (INELgtZERO && col.multNTracksPVeta1() < 1) {
      return;
    }
    if (fabs(col.posZ()) > vtxZsel) {
      return;
    }

    // verify that PV has collision
    if (!col.has_mcCollision())
      return;

    auto mcCollision = col.mcCollision_as<soa::Join<aod::McCollisions, aod::McCollsExtra>>();
    if (mcCollision.hasRecoCollision() < 1)
      return; // total paranoia mode: on

    // Profiles
    if (useZeqInProfiles && do2Dplots) {
      histos.fill(HIST("multiplicityQa/h2dPVsVsFV0"), col.multZeqFV0A(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0"), col.multZeqFT0A() + col.multZeqFT0C(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0A"), col.multZeqFT0A(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0C"), col.multZeqFT0C(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFDD"), col.multZeqFDDA() + col.multZeqFDDC(), mcCollision.hasRecoCollision());
    } else {
      histos.fill(HIST("multiplicityQa/h2dPVsVsFV0"), col.multFV0A(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0"), col.multFT0A() + col.multFT0C(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0A"), col.multFT0A(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFT0C"), col.multFT0C(), mcCollision.hasRecoCollision());
      histos.fill(HIST("multiplicityQa/h2dPVsVsFDD"), col.multFDDA() + col.multFDDC(), mcCollision.hasRecoCollision());
    }
  }
  PROCESS_SWITCH(MultiplicityQa, processCollisions, "per-collision analysis", true);
  PROCESS_SWITCH(MultiplicityQa, processBCs, "per-BC analysis", false);
  PROCESS_SWITCH(MultiplicityQa, processCollisionsPVChecks, "do PV contributors check", false);
  PROCESS_SWITCH(MultiplicityQa, processMCCollisions, "verify mc collisions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiplicityQa>(cfgc)};
}
