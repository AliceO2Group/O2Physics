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

static const std::vector<std::string> mmObjectsNames{"kHtPt", "kHmTrk", "kHmTrkTns", "kHmTrk&kHmTrkTns"};

struct multFilter {
  enum { kLeadingPtTrack = 0,
         kHighTrackMult,
         kHighTrackMultTrans,
         kHighTrackMultOverlap,
         kNtriggersMM };

  // event selection cuts
  Configurable<float> selectionLeadingPtTrack{"selectionLeadingPtTrack", 5., "Minimum track pT leading threshold"};
  Configurable<float> selectionHighTrackMult{"selectionHighTrackMult", 24., "Minimum charged particle multiplicity threshold"};
  Configurable<float> selectionHighTrackMultTrans{"selectionHighTrackMultTrans", 6., "Minimum charged particle multiplicity (in transverse region) threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
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
    multiplicity.add("fTrackMultVsV0A", "charged particle multiplicity", HistType::kTH2F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}, {200, -0.5, +16999.5, "sum AmpFV0"}});
    multiplicity.add("fTrackMultVsT0A", "charged particle multiplicity", HistType::kTH2F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}, {100, -0.5, +299.5, "sum AmpFT0"}});
    multiplicity.add("fTrackMultTrans", "charged particle multiplicity in transverse region", HistType::kTH1F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8, transverse region)"}});
    multiplicity.add("fLeadingTrackPtSelected", "pT of selected high pT tracks", HistType::kTH1F, {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"}});
    multiplicity.add("fTrackMultSelected", "charged particle multiplicity of the selected events", HistType::kTH1F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}});
    multiplicity.add("fTrackMultTransSelected", "charged particle multiplicity  (in the transverse region) of the selected events", HistType::kTH1F, {{200, -0.5, +199.5, "number of tracks (|#eta|<0.8)"}});
    multiplicity.add("fDeltaPhiTrans", "delta phi distribution in transverse region", HistType::kTH1F, {{64, -M_PI / 2.0, +3.0 * M_PI / 2.0, "number of tracks (|#eta|<0.8, trans reg.)"}});

    std::array<std::string, 2> eventTitles = {"all", "rejected"};

    auto scalers{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalers->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalers->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
  }

  // declare filters on tracks and charged jets
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackEtaCut) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::pt > cfgTrackLowPtCut);

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
  float computeDeltaPhi(float phia, float phib,
                        float rangeMin = -M_PI / 2.0, float rangeMax = 3.0 * M_PI / 2.0)
  {
    float dphi = -999;
    if (phia < 0) {
      phia += 2 * M_PI;
    } else if (phia > 2 * M_PI) {
      phia -= 2 * M_PI;
    }
    if (phib < 0) {
      phib += 2 * M_PI;
    } else if (phib > 2 * M_PI) {
      phib -= 2 * M_PI;
    }
    dphi = phib - phia;
    if (dphi < rangeMin) {
      dphi += 2 * M_PI;
    } else if (dphi > rangeMax) {
      dphi -= 2 * M_PI;
    }
    return dphi;
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {

    bool keepEvent[kNtriggersMM]{false};

    multiplicity.fill(HIST("fProcessedEvents"), 0);

    // V0A signal
    float sumAmpFT0 = 0;
    float sumAmpFV0 = 0;
    int innerFV0 = 24;
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        for (auto amplitude : ft0.amplitudeA()) {
          sumAmpFT0 += amplitude;
        }
        for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
          if (int(fv0.channel()[ich]) > innerFV0)
            continue;
          sumAmpFV0 += fv0.amplitude()[ich];
        }
      }
    }

    //
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    int multTrack = 0;
    float flPt = 0; // leading pT
    float flPhi = 0;
    // Check whether there is a leading pT track
    for (auto& track : tracks) { // start loop over tracks
      multTrack++;
      multiplicity.fill(HIST("hdNdeta"), track.eta());
      if (flPt < track.pt()) {
        flPt = track.pt();
        flPhi = track.phi();
      }
    }

    multiplicity.fill(HIST("fLeadingTrackPt"), flPt);
    multiplicity.fill(HIST("fTrackMult"), multTrack);
    multiplicity.fill(HIST("fTrackMultVsV0A"), multTrack, sumAmpFV0);
    multiplicity.fill(HIST("fTrackMultVsT0A"), multTrack, sumAmpFT0);

    // Check whether this event has a leading track candidate
    int multTrackTrans = 0;
    if (flPt >= selectionLeadingPtTrack) {
      // compute Nch transverse
      for (auto& track : tracks) { // start loop over tracks
        float DPhi = computeDeltaPhi(track.phi(), flPhi);
        if (!(TMath::Abs(DPhi) < M_PI / 3.0 || TMath::Abs(DPhi - M_PI) < M_PI / 3.0)) { // transverse side
          multTrackTrans++;
          multiplicity.fill(HIST("fDeltaPhiTrans"), DPhi);
        }
      }

      multiplicity.fill(HIST("fTrackMultTrans"), multTrackTrans);
      if (multTrackTrans >= selectionHighTrackMultTrans) {
        keepEvent[kHighTrackMultTrans] = true; // accepted HM events
        multiplicity.fill(HIST("fTrackMultTransSelected"), multTrackTrans);
        if (keepEvent[kHighTrackMult] && keepEvent[kHighTrackMultTrans]) {
          keepEvent[kHighTrackMultOverlap] = true;
        }
      }
      multiplicity.fill(HIST("fLeadingTrackPtSelected"), flPt); // track pT which passed the cut
      keepEvent[kLeadingPtTrack] = true;
    }

    // Check whether this is a high multiplicity event
    if (multTrack >= selectionHighTrackMult) {
      keepEvent[kHighTrackMult] = true; // accepted HM events
      multiplicity.fill(HIST("fTrackMultSelected"), multTrack);
    }
    tags(keepEvent[kLeadingPtTrack], keepEvent[kHighTrackMult], keepEvent[kHighTrackMultTrans], keepEvent[kHighTrackMultOverlap]);

    if (!keepEvent[kLeadingPtTrack] && !keepEvent[kHighTrackMult] && !keepEvent[kHighTrackMultTrans]) {
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
