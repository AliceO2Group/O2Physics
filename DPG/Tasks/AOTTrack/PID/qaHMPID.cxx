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

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct pidHmpidQa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsP{"nBinsP", 500, "Number of momentum bins"};
  Configurable<float> minP{"minP", 0.01f, "Minimum momentum plotted (GeV/c)"};
  Configurable<float> maxP{"maxP", 10.f, "Maximum momentum plotted (GeV/c)"};
  Configurable<float> maxDCA{"maxDCA", 3.f, "Maximum DCA xy use for the plot (cm)"};
  Configurable<float> maxDistance{"maxDistance", 5.f, "Maximum HMPID distance between the track and the cluster (cm)"};
  Configurable<float> minCharge{"minCharge", 120.f, "Minimum HMPID charge collected in the cluster"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec momAxis{nBinsP, minP, maxP};
    histos.add("qa/signalvsP", ";#it{p} (GeV/#it{c});Cherenkov angle (rad)", kTH2F, {momAxis, {1000, 0, 1}});
    histos.add("distance/selected", ";HMPID distance", kTH1F, {{100, 0, 20}});
    histos.add("distance/nonselected", ";HMPID distance", kTH1F, {{100, 0, 20}});
    histos.add("qmip/selected", ";HMPID mip charge (ADC)", kTH1F, {{100, 0, 4000}});
    histos.add("qmip/nonselected", ";HMPID mip charge (ADC)", kTH1F, {{100, 0, 4000}});
    histos.add("nphotons/selected", ";HMPID number of detected photons", kTH1F, {{100, 0, 1000}});
    histos.add("nphotons/nonselected", ";HMPID number of detected photons", kTH1F, {{100, 0, 1000}});
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  void process(const TrackCandidates& tracks,
               const aod::HMPIDs& hmpids,
               const aod::Collisions& colls)
  {
    for (const auto& t : hmpids) {

      if (t.track_as<TrackCandidates>().isGlobalTrack() != (uint8_t) true) {
        continue;
      }
      if (abs(t.track_as<TrackCandidates>().dcaXY()) > maxDCA) {
        continue;
      }
      histos.fill(HIST("distance/nonselected"), t.hmpidDistance());
      histos.fill(HIST("qmip/nonselected"), t.hmpidQMip());
      histos.fill(HIST("nphotons/nonselected"), t.hmpidNPhotons());
      if (t.hmpidDistance() > maxDistance) {
        continue;
      }
      if (t.hmpidQMip() < minCharge) {
        continue;
      }
      histos.fill(HIST("distance/selected"), t.hmpidDistance());
      histos.fill(HIST("qmip/selected"), t.hmpidQMip());
      histos.fill(HIST("nphotons/selected"), t.hmpidNPhotons());
      histos.fill(HIST("qa/signalvsP"), t.track_as<TrackCandidates>().p(), t.hmpidSignal());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidQa>(cfg)}; }
