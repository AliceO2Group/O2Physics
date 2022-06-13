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
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ALICE3MultTask {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  Configurable<float> MinEta{"MinEta", -0.8f, "Minimum eta in range"};
  Configurable<float> MaxEta{"MaxEta", 0.8f, "Maximum eta in range"};
  Configurable<float> MaxMult{"MaxMult", 1000.f, "Maximum multiplicity in range"};
  Configurable<float> MinContrib{"MinContrib", 2.f, "Minimum number of contributors to the PV"};
  Configurable<float> MaxDCA{"MaxDCA", 25.f, "Max DCAxy and DCAz for counted tracks"};

  void init(InitContext&)
  {
    const AxisSpec axisMult{MaxMult.value > 10000.f ? 10000 : (int)MaxMult, 0, MaxMult, "Reconstructed tracks"};
    TString tit = Form("%.3f < #it{#eta} < %.3f", MinEta.value, MaxEta.value);
    histos.add("multiplicity/numberOfTracksNoCut", tit, kTH1D, {axisMult});
    histos.add("multiplicity/numberOfTracks", tit, kTH1D, {axisMult});
  }

  int nevs = 0;
  void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks)
  {
    int nTracks = 0;
    for (const auto& track : tracks) {
      if (track.eta() < MinEta || track.eta() > MaxEta) {
        continue;
      }
      if (abs(track.dcaXY()) > MaxDCA || abs(track.dcaZ()) > MaxDCA) {
        continue;
      }
      nTracks++;
    }
    LOG(info) << nevs++ << ") Event " << collision.globalIndex() << " has " << nTracks << " tracks";
    histos.fill(HIST("multiplicity/numberOfTracksNoCut"), nTracks);
    if (collision.numContrib() < MinContrib) {
      return;
    }
    histos.fill(HIST("multiplicity/numberOfTracks"), nTracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ALICE3MultTask>(cfgc, TaskName{"alice3-qa-multiplicity"})};
}
