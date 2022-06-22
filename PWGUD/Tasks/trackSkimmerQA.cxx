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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define mmuon 0.1056583745

struct TrackSkimmerQa {
  HistogramRegistry registry{
    "registry",
    {
      {"Timings/TrackTimeRes", ";Track time resolution, ns;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Kine/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
      {"Kine/Eta", ";#eta;", {HistType::kTH1D, {{100, -4., 4.}}}}
    }
  };

  void process(o2::aod::SkimmedMuons const& muonTracks)
  {
    for (const auto& track : muonTracks) {
      TLorentzVector p;
      p.SetXYZM(track.px(), track.py(), track.pz(), mmuon);
      registry.fill(HIST("Timings/TrackTimeRes"), track.trackTimeRes());
      registry.fill(HIST("Kine/Pt"), p.Pt());
      registry.fill(HIST("Kine/Eta"), p.Eta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrackSkimmerQa>(cfgc)};
}
