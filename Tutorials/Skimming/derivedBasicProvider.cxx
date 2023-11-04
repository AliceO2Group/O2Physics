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
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataModel/DerivedExampleTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct DerivedBasicProvider {
  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> minTPCNClsCrossedRows{"minTPCNClsCrossedRows", 80, "min TPC crossed rows"};
  Configurable<float> minPt{"minPt", 2.0, "min pT to save"};
  Configurable<float> maxDCA{"maxDCA", 0.1, "max DCA"};
  Configurable<float> etaWindow{"etaWindow", 0.8, "eta window"};
  Configurable<bool> skipUninterestingEvents{"skipUninterestingEvents", true, "skip collisions without particle of interest"};

  // This marks that this task produces a standard derived table
  Produces<aod::DrCollisions> outputCollisions;
  Produces<aod::DrTracks> outputTracks;

  // Look at primary tracks only
  Filter trackFilter = nabs(aod::track::dcaXY) < maxDCA && nabs(aod::track::eta) < etaWindow && aod::track::pt > minPt;

  // This is an example of a convenient declaration of "using"
  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>; // do not forget this!

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
  }

  void process(aod::Collision const& collision, myFilteredTracks const& tracks)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    if (tracks.size() < 1 && skipUninterestingEvents)
      return;
    bool interestingEvent = false;
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < minTPCNClsCrossedRows)
        continue; // remove badly tracked
      interestingEvent = true;
    }
    if (!interestingEvent && skipUninterestingEvents)
      return;
    outputCollisions(collision.posZ());
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < minTPCNClsCrossedRows)
        continue; // remove badly tracked
      histos.get<TH1>(HIST("ptHistogram"))->Fill(track.pt());
      outputTracks(outputCollisions.lastIndex(), track.pt(), track.eta(), track.phi()); // all that I need for posterior analysis!
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedBasicProvider>(cfgc, TaskName{"derived-basic-provider"})};
  return workflow;
}
