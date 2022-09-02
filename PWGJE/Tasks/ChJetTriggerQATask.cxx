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

// jet Trigger QA Task
//
// Author: Filip Krizek

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "EventFiltering/filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// What this task should do
// Event by event fill
// 1) pT spectrum of tracks in TPC volume
// 2) pT spectrum of jets in fiducial volume
// 3) leading jet pT  versus leading track pT  both in TPC volume
// We want output from
// a) minimum bias events
// b) from events selected by EPN
// It would be good to run it for reveral jet radii  e.g. 0.2, 0.4, 0.6

struct ChJetTriggerQATask {

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0, "Accepted z-vertex range"};
  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9, "Full eta range"};                       //without fiducial cut
  Configurable<float> cfgFiducialVolume{"cfgFiducialVolume", 0.5, "Eta range for charged jets"}; //fiducial cut
  Configurable<int> bTriggerDecision{"bTriggerDecision", 0, "Charged Jet Trigger Decision Seleection"};

  //HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry spectra{
    "spectra",
    {
      {"fTrackPtIncl", "inclusive track pT", {HistType::kTH1F, {{200, 0., +200.}}}},                                      //
      {"fJetChPtIncl", "inclusive charged jet pT", {HistType::kTH1F, {{200, 0., +200.}}}},                                //
      {"fLeadJetChPtVsLeadingTrack", "inclusive charged jet pT", {HistType::kTH2F, {{200, 0., +200.}, {200, 0., +200.}}}} //
    }                                                                                                                     //
  };

  TrackSelection globalTracks;
  void init(o2::framework::InitContext&)
  {

    globalTracks = getGlobalTrackSelection();
    globalTracks.SetEtaRange(-.9, .9);
    //spectra.add("fTrackPtIncl", "inclusive track pT", HistType::kTH1F, {{200, 0., +200., "track #it{p}_{T} (GeV/#it{c})"}});  //in TPC volume
    //spectra.add("fJetChPtIncl", "inclusive charged jet pT", HistType::kTH1F, {{200, 0., +200., "charged jet #it{p}_{T} (GeV/#it{c})"}}); //in fiducial volume
    //spectra.add("fLeadJetChPtVsLeadingTrack", "inclusive charged jet pT", HistType::kTH2F, {{200, 0., +200.} "charged jet #it{p}_{T} (GeV/#it{c})"}, {{200, 0., +200.} "track #it{p}_{T} (GeV/#it{c})"}); //in TPC volume
    //spectra.add("fLeadJetChPtVsLeadingTrack", "inclusive charged jet pT", HistType::kTH2F, {{200, 0., +200.},{200, 0., +200.}}); //in TPC volume
  }

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVertexCut) && (aod::filtering::hasJetChHighPt >= bTriggerDecision);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>>::iterator const& collision, TrackCandidates const& tracks, aod::Jets const& jets)
  {

    //Find leading jet pT in full TPC volume
    double leadingJetPt = -1.0;
    double leadingTrackPt = -1.0;

    for (auto& jet : jets) {
      if (fabs(jet.eta()) < cfgTPCVolume) {
        if (jet.pt() > leadingJetPt) {
          leadingJetPt = jet.pt();
        }
      }
    }

    //Find leading track pT in full TPC volume
    for (auto& trk : tracks) {
      if (!globalTracks.IsSelected(trk))
        continue;
      if (fabs(trk.eta()) < cfgTPCVolume) {
        if (trk.pt() > leadingTrackPt) {
          leadingTrackPt = trk.pt();
        }
      }
    }

    spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt, leadingJetPt); //leading jet pT versus leading track pT

    //--------------------------------------------------------------
    // Inclusive Track pT spectrum in TPC volume
    for (auto& trk : tracks) {
      if (!globalTracks.IsSelected(trk))
        continue;
      if (fabs(trk.eta()) < cfgTPCVolume) {
        spectra.fill(HIST("fTrackPtIncl"), trk.pt());
      }
    }

    // Inclusive Jet pT spectrum in Fiducial volume
    for (auto& jet : jets) {
      if (fabs(jet.eta()) < cfgFiducialVolume) {
        spectra.fill(HIST("fJetChPtIncl"), jet.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{
    adaptAnalysisTask<ChJetTriggerQATask>(cfgc, TaskName{"jet-charged-trigger-qa"})};
}
