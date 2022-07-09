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
//FK #include "Common/DataModel/PIDResponse.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> highPtObjectsNames{"HighPtTrack", "JetChHighPt"};

struct jetFilter {
  enum { kHighPtTrack = 0,
         kJetChHighPt,
         kHighPtObjects };

  //event selection cuts
  Configurable<float> selectionHighPtTrack{"selectionHighPtTrack", 10., "Minimum track pT trigger threshold"};       //we want to keep all events having a track with pT above this
  Configurable<float> selectionJetChHighPt{"selectionJetChHighPt", 40., "Minimum charged jet pT trigger threshold"}; //we want to keep all events having a charged jet with pT above this

  Produces<aod::JetFilters> tags;

  //acceptance cuts
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> cfgTrackEtaCut{"cfgTrackEtaCut", 0.9f, "Eta range for tracks"};
  Configurable<float> cfgTrackLowPtCut{"cfgTrackLowPtCut", 0.1f, "Minimum constituent pT"};

  Configurable<float> cfgJetChEtaCut{"cfgJetChEtaCut", 0.5f, "Eta range for charged jets"};

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    spectra.add("fTrackPtSelected", "pT of selected high pT tracks", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    spectra.add("fJetChPtSelected", "pT of selected high pT charged jets", HistType::kTH1F, {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"}});

    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add("fProcessedEvents", ";;Number of filtered events", HistType::kTH1F, {{kHighPtObjects, -0.5, kHighPtObjects - 0.5}}))};
    for (uint32_t iS{1}; iS <= highPtObjectsNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, highPtObjectsNames[iS - 1].data());
    }
  }

  //declare filters on tracks and charged jets
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > cfgTrackLowPtCut);
  Filter jetChFilter = (nabs(aod::jet::eta) < cfgJetChEtaCut);

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks, aod::Jets const& jets)
  {
    // collision process loop
    bool keepEvent[kHighPtObjects]{false};
    //
    spectra.fill(HIST("fCollZpos"), collision.posZ());

    //Check whether there is a high pT track
    for (auto& track : tracks) { // start loop over tracks
      if (track.pt() >= selectionHighPtTrack) {
        spectra.fill(HIST("fTrackPtSelected"), track.pt()); //track pT which passed the cut
        keepEvent[kHighPtTrack] = true;
        break;
      }
    }

    //Check whether there is a high pT charged jet
    for (auto& jet : jets) { // start loop over charged jets
      if (jet.pt() >= selectionJetChHighPt) {
        spectra.fill(HIST("fJetChPtSelected"), jet.pt()); //charged jet pT
        keepEvent[kJetChHighPt] = true;
        break;
      }
    }

    //count events which passed the selections
    for (int iDecision{0}; iDecision < kHighPtObjects; ++iDecision) {
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }
    tags(keepEvent[0], keepEvent[1]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  return WorkflowSpec{
    adaptAnalysisTask<jetFilter>(cfg)};
}
