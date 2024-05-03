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
  Configurable<bool> enabledEdxPerID{"enabledEdxPerID", false, "Flag to produce dE/dx per particle ID histograms"};
  Configurable<int> trdSelection{"trdSelection", 0, "Flag for the TRD selection: -1 no TRD, 0 no selection, 1 TRD"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxEta{"maxEta", 0.8, "Maximum eta in range"};
  Configurable<float> maxITSChi2{"maxITSChi2", 36, "Maximum chi2 in range"};
  Configurable<float> maxTPCChi2{"maxTPCChi2", 4, "Maximum chi2 in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<int> minITSCls{"minITSCls", 6, "Minimum number of ITS clusters"};
  Configurable<int> pidInTracking{"pidInTracking", -1, "PID in tracking"};
  Configurable<int> minTPCCls{"minTPCCls", 0, "Minimum number of TPC clusters"};
  Configurable<int> minTPCClsFindable{"minTPCClsFindable", 0, "Minimum number of TPC findable clusters"};
  Configurable<float> minCrossedRowsOverFindableCls{"minCrossedRowsOverFindableCls", 0.8f, "Minimum number of TPC found/findable clusters"};
  ConfigurableAxis dEdxBins{"dEdxBins", {5000, 0.f, 5000.f}, "Binning in dE/dx"};
  Configurable<float> minTPCNcls{"minTPCNcls", 0.f, "Minimum number or TPC Clusters for tracks"};
  Configurable<float> minNClsCrossedRows{"minNClsCrossedRows", 70.f, "Minimum number or TPC crossed rows for tracks"};
  Configurable<bool> enableITS{"enableITS", false, "Enable ITS plots"};

  std::array<std::shared_ptr<TH2>, 2> hdedx;
  std::array<std::array<std::shared_ptr<TH2>, 9>, 2> hdedx_perID;
  std::array<std::shared_ptr<TH3>, 2> hdedx_its;
  std::array<std::array<std::shared_ptr<TH3>, 9>, 2> hdedx_perID_its;
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

    // Event properties
    auto h = histos.add<TH1>("event/evsel", "", kTH1D, {{10, 0.5, 10.5, "Ev. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Passed ev. sel.");
    h->GetXaxis()->SetBinLabel(3, "Passed vtx Z");

    // Track properties
    h = histos.add<TH1>("trksel", "", kTH1D, {{10, 0.5, 10.5, "Trk. Sel."}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, Form("Has > %i ITS clusters", minITSCls.value));
    h->GetXaxis()->SetBinLabel(3, Form("Has > %i TPC clusters findable", minTPCClsFindable.value));
    h->GetXaxis()->SetBinLabel(4, Form("Has > %f TPC clusters found", minTPCNcls.value));
    h->GetXaxis()->SetBinLabel(5, Form("Has > %f Found/Findable Ratio", minCrossedRowsOverFindableCls.value));
    h->GetXaxis()->SetBinLabel(6, Form("Has > %f Xrows", minNClsCrossedRows.value));
    h->GetXaxis()->SetBinLabel(7, "All PID in trk");
    if (pidInTracking == -1) {
      h->GetXaxis()->SetBinLabel(7, Form("PID in trk %i", pidInTracking.value));
    }
    h->GetXaxis()->SetBinLabel(8, "no TRD sel.");
    if (trdSelection == -1) {
      h->GetXaxis()->SetBinLabel(8, "has no TRD");
    } else if (trdSelection == 1) {
      h->GetXaxis()->SetBinLabel(8, "has TRD");
    }

    histos.add("event/vertexz", "", kTH1D, {vtxZAxis});

    const AxisSpec itsAxis{20, -0.5, 19.5, "<ITS cluster size>"};
    if (enableITS) {
      if (enabledEdxPerID) {
        for (int i = 0; i < 9; i++) {
          hdedx_perID_its[0][i] = histos.add<TH3>(Form("event/tpcsignalPos_%i", i), "", kTH3F, {pAxis, dedxAxis, itsAxis});
          hdedx_perID_its[1][i] = histos.add<TH3>(Form("event/tpcsignalNeg_%i", i), "", kTH3F, {pAxis, dedxAxis, itsAxis});
        }
      } else {
        hdedx_its[0] = histos.add<TH3>("event/tpcsignalPos", "", kTH3F, {pAxis, dedxAxis, itsAxis});
        hdedx_its[1] = histos.add<TH3>("event/tpcsignalNeg", "", kTH3F, {pAxis, dedxAxis, itsAxis});
      }
    } else {
      if (enabledEdxPerID) {
        for (int i = 0; i < 9; i++) {
          hdedx_perID[0][i] = histos.add<TH2>(Form("event/tpcsignalPos_%i", i), "", kTH2D, {pAxis, dedxAxis});
          hdedx_perID[1][i] = histos.add<TH2>(Form("event/tpcsignalNeg_%i", i), "", kTH2D, {pAxis, dedxAxis});
        }
      } else {
        hdedx[0] = histos.add<TH2>("event/tpcsignalPos", "", kTH2D, {pAxis, dedxAxis});
        hdedx[1] = histos.add<TH2>("event/tpcsignalNeg", "", kTH2D, {pAxis, dedxAxis});
      }
    }

    LOG(info) << "QA PID TPC histograms:";
    histos.print();
  }

  template <typename T>
  bool isTrackSelected(const T& track)
  {
    histos.fill(HIST("trksel"), 1);
    if (track.itsNCls() < minITSCls) {
      return false;
    }
    histos.fill(HIST("trksel"), 2);
    if (track.tpcNClsFindable() <= minTPCClsFindable) {
      return false;
    }
    histos.fill(HIST("trksel"), 3);
    if (track.tpcNClsFound() < minTPCNcls) {
      return false;
    }
    histos.fill(HIST("trksel"), 4);
    if (track.tpcCrossedRowsOverFindableCls() < minCrossedRowsOverFindableCls) {
      return false;
    }
    histos.fill(HIST("trksel"), 5);
    if (track.tpcNClsCrossedRows() < minNClsCrossedRows) {
      return false;
    }
    histos.fill(HIST("trksel"), 6);
    if (pidInTracking != -1 && (track.pidForTracking() != std::abs(pidInTracking))) {
      return false;
    }
    histos.fill(HIST("trksel"), 7);
    switch (trdSelection) {
      case 0:
        break;
      case -1:
        if (track.hasTRD()) {
          return false;
        }
        break;
      case 1:
        if (!track.hasTRD()) {
          return false;
        }
        break;
      default:
        LOG(fatal) << "Invalid TRD selection";
    }
    histos.fill(HIST("trksel"), 8);
    return true;
  }

  Filter trackFilterEta = (nabs(aod::track::eta) < maxEta);
  Filter trackFilterITS = (aod::track::itsChi2NCl < maxITSChi2);
  Filter trackFilterTPC = ((aod::track::tpcChi2NCl < maxTPCChi2));
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra>;
  void processEvSel(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                    soa::Filtered<TrackCandidates> const& tracks,
                    aod::BCs const&)
  {
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
      if (!isTrackSelected(t)) {
        continue;
      }
      switch (trdSelection) {
        case 0:
          break;
        case -1:
          if (t.hasTRD()) {
            continue;
          }
          break;
        case 1:
          if (!t.hasTRD()) {
            continue;
          }
          break;
        default:
          LOG(fatal) << "Invalid TRD selection";
      }
      const int signIdx = t.sign() > 0 ? 0 : 1;
      hdedx[signIdx]->Fill(t.tpcInnerParam(), t.tpcSignal());
      if (enabledEdxPerID) {
        hdedx_perID[signIdx][t.pidForTracking()]->Fill(t.tpcInnerParam(), t.tpcSignal());
      }
    }
  }
  PROCESS_SWITCH(tpcPidQaSignal, processEvSel, "Process with event selection", false);

  void processNoEvSel(soa::Filtered<TrackCandidates> const& tracks, aod::Collisions const& collisions)
  {
    histos.fill(HIST("event/evsel"), 1, collisions.size());

    for (const auto& t : tracks) {
      if (!t.has_collision()) {
        continue;
      }
      if (abs(t.collision().posZ()) > 10.f) {
        continue;
      }
      if (!isTrackSelected(t)) {
        continue;
      }

      const int signIdx = t.sign() > 0 ? 0 : 1;
      if (enableITS) {
        float average = 0;
        int n = 0;
        for (int i = 0; i < 7; i++) {
          if (t.itsClsSizeInLayer(i) <= 0) {
            continue;
          }
          average += t.itsClsSizeInLayer(i);
          n++;
        }
        average = n > 0 ? average / n : 0;
        average *= cos(atan(t.tgl()));
        if (enabledEdxPerID) {
          hdedx_perID_its[signIdx][t.pidForTracking()]->Fill(t.tpcInnerParam(), t.tpcSignal(), average);
        } else {
          hdedx_its[signIdx]->Fill(t.tpcInnerParam(), t.tpcSignal(), average);
        }
      } else {
        if (enabledEdxPerID) {
          hdedx_perID[signIdx][t.pidForTracking()]->Fill(t.tpcInnerParam(), t.tpcSignal());
        } else {
          hdedx[signIdx]->Fill(t.tpcInnerParam(), t.tpcSignal());
        }
      }
    }
  }
  PROCESS_SWITCH(tpcPidQaSignal, processNoEvSel, "Process without event selection", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tpcPidQaSignal>(cfgc)}; }
