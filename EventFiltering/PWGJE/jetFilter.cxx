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

// Author: Filip Krizek

#include <TMath.h>
#include <cmath>
#include <string>

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> highPtObjectsNames{"JetChHighPt"};

struct jetFilter {
  enum { kJetChHighPt = 0,
         kHighPtObjects };

  // event selection cuts
  Configurable<float> selectionJetChHighPt{
    "selectionJetChHighPt", 33.,
    "Minimum charged jet pT trigger threshold"}; // we want to keep all events
                                                 // having a charged jet with
                                                 // pT above this

  Produces<aod::JetFilters> tags;

  // acceptance cuts
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0,
                                   "Accepted z-vertex range"};

  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9,
                                   "Full eta range"}; // eta range of TPC
  Configurable<float> cfgJetR{"cfgJetR", 0.6,
                              "jet resolution parameter"}; // jet cone radius
  Configurable<float> cfgJetPtMin{
    "cfgJetPtMin", 0.05,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT

  HistogramRegistry spectra{
    "spectra",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(o2::framework::InitContext&)
  {

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F,
                {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});
    spectra.add(
      "ptphiGoodTracks", "ptphiGoodTracks", HistType::kTH2F,
      {{100, 0., +100., "#it{p}_{T} (GeV)"}, {60, 0, TMath::TwoPi()}});
    spectra.add(
      "ptphiRejectedTracks", "ptphiRejectedTracks", HistType::kTH2F,
      {{100, 0., +100., "#it{p}_{T} (GeV)"}, {60, 0, TMath::TwoPi()}});
    spectra.add("ptetaGoodTracks", "ptetaGoodTracks", HistType::kTH2F,
                {{100, 0., +100., "#it{p}_{T} (GeV)"}, {36, -0.9, 0.9}});
    spectra.add("ptetaRejectedTracks", "ptetaRejectedTracks", HistType::kTH2F,
                {{100, 0., +100., "#it{p}_{T} (GeV)"}, {36, -0.9, 0.9}});

    spectra.add("ptphiJetChSelected",
                "pT of selected high pT charged jets vs phi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi()}});
    spectra.add("ptetaJetChSelected",
                "pT of selected high pT charged jets vs eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {36, -0.9, 0.9}});

    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add(
      "fProcessedEvents", ";;Number of filtered events", HistType::kTH1F,
      {{kHighPtObjects, -0.5, kHighPtObjects - 0.5}}))};
    for (uint32_t iS{1}; iS <= highPtObjectsNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, highPtObjectsNames[iS - 1].data());
    }

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = cfgJetR;
    jetReclusterer.jetEtaMin = -2 * cfgTPCVolume;
    jetReclusterer.jetEtaMax = 2 * cfgTPCVolume;
  }

  // declare filters on tracks
  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgVertexCut;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,
                                    aod::TracksDCA, aod::TrackSelection>;

  void
    process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
            TrackCandidates const& tracks)
  {
    // collision process loop
    bool keepEvent[kHighPtObjects]{false};
    //
    spectra.fill(HIST("fCollZpos"), collision.posZ());

    jetConstituents.clear();
    jetReclustered.clear();

    for (auto& trk : tracks) {
      if (!trk.isQualityTrack()) {
        spectra.fill(HIST("ptphiRejectedTracks"), trk.pt(), trk.phi());
        spectra.fill(HIST("ptetaRejectedTracks"), trk.pt(), trk.eta());
        continue;
      }
      if (fabs(trk.eta()) < cfgTPCVolume) {
        spectra.fill(HIST("ptphiGoodTracks"), trk.pt(), trk.phi());
        spectra.fill(HIST("ptetaGoodTracks"), trk.pt(), trk.eta());

        if (trk.pt() > cfgJetPtMin) { // jet constituents
          fillConstituents(trk,
                           jetConstituents); // ./PWGJE/Core/JetFinder.h
                                             // recombination scheme is assumed
                                             // to be Escheme with pion mass
        }
      }
    }

    // Reconstruct jet from tracks
    fastjet::ClusterSequenceArea clusterSeq(
      jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);

    // Check whether there is a high pT charged jet
    for (auto& jet : jetReclustered) { // start loop over charged jets
      if (fabs(jet.eta()) < 2 * cfgTPCVolume) {
        if (jet.perp() >= selectionJetChHighPt) {
          spectra.fill(HIST("ptphiJetChSelected"), jet.perp(),
                       jet.phi()); // charged jet pT vs phi
          spectra.fill(HIST("ptetaJetChSelected"), jet.perp(),
                       jet.eta()); // charged jet pT vs eta
          keepEvent[kJetChHighPt] = true;
          break;
        }
      }
    }

    // count events which passed the selections
    if (fabs(collision.posZ()) > cfgVertexCut)
      keepEvent[kJetChHighPt] = false;

    for (int iDecision{0}; iDecision < kHighPtObjects; ++iDecision) {
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }
    tags(collision, keepEvent[kJetChHighPt]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  return WorkflowSpec{adaptAnalysisTask<jetFilter>(cfg)};
}
