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
// \file multFilter.cxx
// \brief task for selection of high multiplicity events
//
// \author Antonio Ortiz <antonio.ortiz.velasquez@cern.ch>, ICN-UNAM
//
// O2 includes

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> mmObjectsNames{"LeadingPtTrack", "HighTrackMult"};

struct multFilter {
  enum { kLeadingPtTrack = 0,
         kHighTrackMult,
         kNtriggersMM };

  //event selection cuts
  Configurable<float> selectionLeadingPtTrack{"selectionLeadingPtTrack", 8., "Minimum track pT leading threshold"};
  Configurable<float> selectionHighTrackMult{"selectionHighTrackMult", 58., "Minimum charged particle multiplicity threshold"};

  Produces<aod::MultFilters> tags;

  //acceptance cuts
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> cfgTrackEtaCut{"cfgTrackEtaCut", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgTrackLowPtCut{"cfgTrackLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  void init(o2::framework::InitContext&)
  {

    multiplicity.add("fCollZpos", "collision z position", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    multiplicity.add("hdNdeta", "dNdeta", HistType::kTH1F, {{50, -2.5, 2.5, " "}});
    multiplicity.add("fLeadingTrackPt", "pT of high pT tracks", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fTrackMult", "charged particle multiplicity", HistType::kTH1F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}});

    multiplicity.add("fLeadingTrackPtSelected", "pT of selected high pT tracks", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fTrackMultSelected", "charged particle multiplicity of the selected events", HistType::kTH1F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}});

    std::array<std::string, 2> eventTitles = {"all", "rejected"};

    auto scalers{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalers->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalers->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
  }

  //declare filters on tracks and charged jets
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrackLowPtCut);

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks)
  {

    bool keepEvent[kNtriggersMM]{false};

    multiplicity.fill(HIST("fProcessedEvents"), 0);

    //
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    Int_t multTrack = 0;
    Double_t flPt = 0; // leading pT
    //Check whether there is a leading pT track
    for (auto& track : tracks) { // start loop over tracks
      multTrack++;
      multiplicity.fill(HIST("hdNdeta"), track.eta());
      if (flPt < track.pt()) {
        flPt = track.pt();
      }
    }

    multiplicity.fill(HIST("fLeadingTrackPt"), flPt);
    multiplicity.fill(HIST("fTrackMult"), multTrack);
    // Check whether this event has a leading track candidate
    if (flPt >= selectionLeadingPtTrack) {
      multiplicity.fill(HIST("fLeadingTrackPtSelected"), flPt); //track pT which passed the cut
      keepEvent[kLeadingPtTrack] = true;
    }

    //Check whether this is a high multiplicity event
    if (multTrack >= selectionHighTrackMult) {
      keepEvent[kHighTrackMult] = true; // accepted HM events
      multiplicity.fill(HIST("fTrackMultSelected"), multTrack);
    }
    tags(keepEvent[kLeadingPtTrack], keepEvent[kHighTrackMult]);

    if (!keepEvent[kLeadingPtTrack] && !keepEvent[kHighTrackMult]) {
      multiplicity.fill(HIST("fProcessedEvents"), 1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
        if (keepEvent[iTrigger]) {
          multiplicity.fill(HIST("fProcessedEvents"), iTrigger + 2);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  return WorkflowSpec{
    adaptAnalysisTask<multFilter>(cfg)};
}
