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
/// \file   qaPIDTPCSignal.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  24/03/2023
/// \brief  Implementation for QA tasks of the TPC signal (lite task, refer to qaPIDTPC.cxx for full TPC QA for PID)
///

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Task to produce the TPC QA plots
struct tpcPidQaSignal {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> logAxis{"logAxis", 1, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<bool> runPerRunOutput{"runPerRunOutput", false, "Flag to produce run by run output for e.g. lite calibration purposes"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 0.1f, "Downsampling factor for the events for derived data"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  ConfigurableAxis dEdxBins{"dEdxBins", {5000, 0.f, 5000.f}, "Binning in dE/dx"};
  Configurable<float> minTPCNcls{"minTPCNcls", 0.f, "Minimum number or TPC Clusters for tracks"};
  Configurable<float> minNClsCrossedRows{"minNClsCrossedRows", 70.f, "Minimum number or TPC crossed rows for tracks"};

  std::shared_ptr<TH3> hdedx;
  int lastRun = -1;
  unsigned int randomSeed = 0;
  void init(o2::framework::InitContext&)
  {
    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p}/|Z| (GeV/#it{c})"};
    if (logAxis) {
      pAxis.makeLogarithmic();
    }
    const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1D, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");

    histos.add("event/vertexz", "", kTH1D, {vtxZAxis});
    hdedx = histos.add<TH3>("event/tpcsignal", "", kTH3D, {pAxis, dedxAxis, chargeAxis});
    LOG(info) << "QA PID TPC histograms:";
    histos.print();
  }

  Filter trackFilterEta = (nabs(aod::track::eta) < 0.8f);
  Filter trackFilterITS = (aod::track::itsChi2NCl < 36.f);
  Filter trackFilterTPC = ((aod::track::tpcChi2NCl < 4.f));
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra>;
  void processEvSel(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                    soa::Filtered<TrackCandidates> const& tracks,
                    aod::BCs const&)
  {
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return;
    }

    histos.fill(HIST("event/evsel"), 1);
    if (!collision.sel8()) {
      return;
    }

    histos.fill(HIST("event/evsel"), 2);

    if (abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("event/evsel"), 3);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& t : tracks) {
      if (t.itsNCls() < 6) {
        continue;
      }
      if (t.tpcNClsFindable() <= 0) {
        continue;
      }
      if (t.tpcNClsFound() < minTPCNcls) {
        continue;
      }
      if (t.tpcCrossedRowsOverFindableCls() < 0.8f) {
        continue;
      }
      if (t.tpcNClsCrossedRows() < minNClsCrossedRows) {
        continue;
      }

      if (runPerRunOutput.value == true && lastRun != collision.bc().runNumber()) {
        lastRun = collision.bc().runNumber();
        AxisSpec pAxis{nBinsP, minP, maxP, "#it{p}/|Z| (GeV/#it{c})"};
        if (logAxis) {
          pAxis.makeLogarithmic();
        }
        const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
        const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};
        hdedx = histos.add<TH3>(Form("Run%i/tpcsignal", lastRun), "", kTH3D, {pAxis, dedxAxis, chargeAxis});
      }
      hdedx->Fill(t.tpcInnerParam(), t.tpcSignal(), t.sign());
    }
  }
  PROCESS_SWITCH(tpcPidQaSignal, processEvSel, "Process with event selection", false);

  void processNoEvSel(aod::Collision const& collision,
                      soa::Filtered<TrackCandidates> const& tracks)
  {
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return;
    }

    histos.fill(HIST("event/evsel"), 1);
    histos.fill(HIST("event/evsel"), 2);

    if (abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("event/evsel"), 3);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& t : tracks) {
      if (t.itsNCls() < 6) {
        continue;
      }
      if (t.tpcNClsFindable() <= 0) {
        continue;
      }
      if (t.tpcNClsFound() < minTPCNcls) {
        continue;
      }
      if (t.tpcCrossedRowsOverFindableCls() < 0.8f) {
        continue;
      }
      if (t.tpcNClsCrossedRows() < minNClsCrossedRows) {
        continue;
      }
      histos.fill(HIST("event/tpcsignal"), t.tpcInnerParam(), t.tpcSignal(), t.sign());
    }
  }
  PROCESS_SWITCH(tpcPidQaSignal, processNoEvSel, "Process without event selection", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tpcPidQaSignal>(cfgc)}; }
