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
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


struct jetFilter {
  enum { kJetChLowPt = 0,
         kJetChHighPt = 1,
         kTrackPt = 2,
         kTriggerObjects = 3
         };

  enum { kBinAllEvents=0,
	 kBinEventsVz=1,
	 kBinJetChLowPt = 2,
         kBinJetChHighPt = 3,
         kBinTrackPt = 4,
         kBins = 5 
         };


  Produces<aod::JetFilters> tags;

  Configurable<float> cfgJetR{"cfgJetR", 0.6,
                              "trigger jet resolution parameter"}; // jet cone radius

  Configurable<float> trackPtLowThreshold{"trackPtThreshold", 25.0, "low pT threshold for track trigger"};
  Configurable<float> jetPtLowThreshold{"jetPtLowThreshold", 30.0, "threshold for low pT jet trigger"};
  Configurable<float> jetPtHighThreshold{"jetPtHighThreshold", 50.0, "threshold for high pT jet trigger"};

  Configurable<float> cfgZvtx{"cfgZvtx", 10.,
                              "z vertex cut"}; // z vertex cut for full stat. inclusing jet spectra

  Configurable<float> cfgEtaTPC{"cfgEtaTPC", 0.9,
                                "pseudorapidity coverage of TPC"}; // for full stat. inclusive jet spectra

  // jet resolution parameters for inclusive jets in fiducial volume
  ConfigurableAxis cfgJetRadii{"cfgJetRadii", {VARIABLE_WIDTH, 0.2, 0.4, 0.6, 1.}, "jet resolution parameters (KEEP THE LAST PARAMETER AT 1.)"};

  std::vector<double> jetRFidVolume; // pseudorapidity limit for given a jet with given R
  std::vector<int> jetIntR;          // jet radius * 100
  int triggerJetR;

  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", kBins, 0.0f, static_cast<float>(kBins))};

  HistogramRegistry spectra{
    "spectra",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  Filter trackFilter = (nabs(aod::jtrack::eta) < static_cast<float>(cfgEtaTPC)) && (aod::jtrack::pt > 0.15f);

  void init(o2::framework::InitContext&)
  {
    triggerJetR = TMath::Nint(cfgJetR * 100.0f);

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F,
                {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    spectra.add("ptphiJetChSelected_lowptjettrigger", "#it{p}_{T} of selected low pT charged jet trigger vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptphiJetChSelected_highptjettrigger", "#it{p}_{T} of selected high pT charged jet trigger vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptphiTrackSelected_trackpttrigger", "#it{p}_{T} of selected high pT tracks vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptetaJetChSelected_lowptjettrigger", "#it{p}_{T} of selected low pT charged jet trigger vs #eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});

    spectra.add("ptetaJetChSelected_highptjettrigger", "#it{p}_{T} of selected high pT charged jet trigger vs #eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});

    spectra.add("ptetaTrackSelected_trackpttrigger", "#it{p}_{T} of selected high pT tracks vs #eta", HistType::kTH2F,
                {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});



    AxisSpec jetRadiiAxis = {cfgJetRadii, "AKT jet resolution parameters"};
    const AxisSpec axisTrackPt{200, 0., +200., "#it{p}_{T,track} (GeV/#it{c})"};
    const AxisSpec axisJetPt{1020, -20., +1000., "#it{p}_{T,jet} (GeV/#it{c})"};
    const AxisSpec axisEta{100, -1., +1., "#eta"};
    const AxisSpec axisPhi{100, 0., TMath::TwoPi(), "#varphi"};

    // inclusive jet spectra in TPC fiducial volume for events with zvtx < 10 cm
    for (unsigned int ir = 0; ir < jetRadiiAxis.binEdges.size() - 1; ir++) {
      jetRFidVolume.push_back(cfgEtaTPC - jetRadiiAxis.binEdges[ir]);
      jetIntR.push_back(TMath::Nint(100 * jetRadiiAxis.binEdges[ir]));
    }
    spectra.add("hEtaVsPtTracksInclusive", "#eta of tracks |#eta| < 0.9;",
                {HistType::kTH2F, {axisTrackPt, axisEta}});

    spectra.add("hPhiVsPtTracksInclusive", "#varphi of tracks |#eta| < 0.9;",
                {HistType::kTH2F, {axisTrackPt, axisPhi}});

    spectra.add("hPtAKTJetsInclusive", "#it{p}_{T} of AKT charged jets in |#eta| < 0.9 - #it{R};",
                {HistType::kTH2F, {jetRadiiAxis, axisJetPt}});

    spectra.add("hPtAKTJetsInclusiveBgSubtr", "#it{p}_{T} of AKT charged jets in |#eta| < 0.9 - #it{R};",
                {HistType::kTH2F, {jetRadiiAxis, axisJetPt}});

    spectra.add("hEtaAKTJetsInclusive", "#eta of AKT charged jets with #it{p}_{T} > 10 GeV in |#eta| < 0.9 - #it{R};",
                {HistType::kTH2F, {jetRadiiAxis, axisEta}});

    spectra.add("hPhiAKTJetsInclusive", "#varphi of AKT charged jets with #it{p}_{T} > 10 GeV in |#eta| < 0.9 - #it{R};",
                {HistType::kTH2F, {jetRadiiAxis, axisPhi}});

    spectra.add("hRho", "Underlying event density #rho", HistType::kTH1F,
                {{200, 0., +20., "#rho (GeV/#it{c})"}});

    hProcessedEvents->GetXaxis()->SetBinLabel(kBinAllEvents+1, "Processed events");
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinEventsVz+1, "Processed events with |vz| < 10 cm");
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinJetChLowPt+1, "JetChLowPt");
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinJetChHighPt+1, "JetChHighPt");
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinTrackPt+1, "TrackPt");
  }

  // declare filters on tracks
  // Filter collisionFilter = nabs(aod::jcollision::posZ) < cfgVertexCut;

  // Filter jetRadiusSelection = o2::aod::jet::r == nround(cfgJetR.node() * 100.0f);
  // FK//using filteredJets = o2::soa::Filtered<o2::aod::ChargedJets>;

  // void process(aod::JCollision const& collision, o2::aod::ChargedJets const& jets)

  template <bool withRho, typename T, typename U, typename D>
  void doTriggering(T const& collision, U const& jets, D const& tracks)
  {

    // collision process loop
    bool keepEvent[kTriggerObjects]{false};
    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt], keepEvent[kTrackPt]);
      return;
    }

    spectra.fill(HIST("fCollZpos"), collision.posZ());
    hProcessedEvents->Fill(static_cast<float>(kBinAllEvents) + 1e-5); // all minimum bias events

    // FILL SPECTRA OF INCLUSIVE JETS IN FIDUCIAL VOLUME
    if (TMath::Abs(collision.posZ()) < cfgZvtx) {
      hProcessedEvents->Fill(static_cast<float>(kBinEventsVz) + 1e-5); // minimum bias events |z_vtx|<10 cm
      if constexpr (withRho) {
        spectra.fill(HIST("hRho"), collision.rho());
      }

      for (const auto& jet : jets) { // jets are ordered by pT
        for (unsigned int ir = 0; ir < jetIntR.size(); ir++) {
          if (jet.r() == jetIntR[ir]) {
            if (TMath::Abs(jet.eta()) < jetRFidVolume[ir]) {
              float jetr = (jet.r() / 100. + 1e-5);
              spectra.fill(HIST("hPtAKTJetsInclusive"), jetr, jet.pt());
              if constexpr (withRho) {
                spectra.fill(HIST("hPtAKTJetsInclusiveBgSubtr"), jetr, jet.pt() - (collision.rho() * jet.area()));
              }

              if (jet.pt() > 10.) {
                spectra.fill(HIST("hEtaAKTJetsInclusive"), jetr, jet.eta());
                spectra.fill(HIST("hPhiAKTJetsInclusive"), jetr, jet.phi());
              }
            }
          }
        }
      }
    }

    for (const auto& jet : jets) { // jets are ordered by pT
      if (jet.r() != triggerJetR)
        continue;

      if (jet.pt() >= jetPtLowThreshold) {
        spectra.fill(HIST("ptphiJetChSelected_lowptjettrigger"), jet.pt(), jet.phi()); // charged jet pT vs phi
        spectra.fill(HIST("ptetaJetChSelected_lowptjettrigger"), jet.pt(), jet.eta()); // charged jet pT vs eta
        keepEvent[kJetChLowPt] = true;
      }
      if (jet.pt() >= jetPtHighThreshold) {
        spectra.fill(HIST("ptphiJetChSelected_highptjettrigger"), jet.pt(), jet.phi()); // charged jet pT vs phi
        spectra.fill(HIST("ptetaJetChSelected_highptjettrigger"), jet.pt(), jet.eta()); // charged jet pT vs eta
        keepEvent[kJetChHighPt] = true;
      }
      break; // only looks at the highest pT jet in the event
    }

    float leadingTrackPt  = -1.;
    float leadingTrackPhi = -100.;
    float leadingTrackEta = -100.;
    for (const auto& track : tracks) { //search for the leading track
      spectra.fill(HIST("hPhiVsPtTracksInclusive"), track.pt(), track.phi());
      spectra.fill(HIST("hEtaVsPtTracksInclusive"), track.pt(), track.eta());
      if(track.pt() > leadingTrackPt){
	leadingTrackPt = track.pt();
        leadingTrackPhi = track.phi();
        leadingTrackEta = track.eta();
      }
    }
    if(leadingTrackPt > trackPtLowThreshold){
      keepEvent[kTrackPt] = true;
      spectra.fill(HIST("ptphiTrackSelected_trackpttrigger"),leadingTrackPt,leadingTrackPhi);
      spectra.fill(HIST("ptetaTrackSelected_trackpttrigger"),leadingTrackPt,leadingTrackEta);
    }

    if(keepEvent[kJetChLowPt]) hProcessedEvents->Fill(static_cast<float>(kBinJetChLowPt) + 1e-5);
    if(keepEvent[kJetChHighPt]) hProcessedEvents->Fill(static_cast<float>(kBinJetChHighPt) + 1e-5);
    if(keepEvent[kTrackPt]) hProcessedEvents->Fill(static_cast<float>(kBinTrackPt) + 1e-5);

    tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt], keepEvent[kTrackPt]);
  }

  void processWithoutRho(soa::Join<JetCollisions, aod::EvSels>::iterator const& collision, o2::aod::ChargedJets const& jets, soa::Filtered<JetTracks> const& tracks) 
  {
    doTriggering<false>(collision, jets, tracks);
  }
  PROCESS_SWITCH(jetFilter, processWithoutRho, "Do charged jet triggering without background estimation for filling histograms", true);

  void processWithRho(soa::Join<JetCollisions, aod::BkgChargedRhos, aod::EvSels>::iterator const& collision, o2::aod::ChargedJets const& jets, soa::Filtered<JetTracks> const& tracks) 
  {
    doTriggering<true>(collision, jets, tracks);
  }
  PROCESS_SWITCH(jetFilter, processWithRho, "Do charged jet triggering with background estimation for filling histograms", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<jetFilter>(cfg)};
}
