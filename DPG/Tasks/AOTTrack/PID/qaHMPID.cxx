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
    histos.add("hmpidSignal", "hmpidSignal", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidXTrack", "hmpidXTrack", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidYTrack", "hmpidYTrack", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidXMip", "hmpidXMip", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidYMip", "hmpidYMip", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidNPhotons", "hmpidNPhotons", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidQMip", "hmpidQMip", kTH1F, {{300, 200, 3200}});
    histos.add("hmpidClusSize", "hmpidClusSize", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidMom", "hmpidMom", kTH1F, {{300, 0, 3000}});
    histos.add("hmpidPhotsCharge", "hmpidPhotsCharge", kTH2F, {{300, 0, 3000}, {300, 0, 3000}});
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  void process(const aod::HMPIDs& hmpids,
               const TrackCandidates& tracks,
               const aod::Collisions& colls)
  {
    for (const auto& t : hmpids) {
      if (t.track_as<TrackCandidates>().isGlobalTrack() != (uint8_t) true) {
        continue;
      }
      if (abs(t.track_as<TrackCandidates>().dcaXY()) > maxDCA) {
        continue;
      }

      histos.fill(HIST("hmpidSignal"), t.hmpidSignal());
      histos.fill(HIST("hmpidXTrack"), t.hmpidXTrack());
      histos.fill(HIST("hmpidYTrack"), t.hmpidYTrack());
      histos.fill(HIST("hmpidXMip"), t.hmpidXMip());
      histos.fill(HIST("hmpidYMip"), t.hmpidYMip());
      histos.fill(HIST("hmpidNPhotons"), t.hmpidNPhotons());
      histos.fill(HIST("hmpidQMip"), t.hmpidQMip());
      histos.fill(HIST("hmpidClusSize"), t.hmpidClusSize());
      histos.fill(HIST("hmpidMom"), t.hmpidMom());
      for (int i = 0; i < 10; i++) {
        histos.fill(HIST("hmpidPhotsCharge"), t.hmpidMom(), i);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidQa>(cfg)}; }
