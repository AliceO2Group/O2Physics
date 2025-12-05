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

#include "../filterTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TMath.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct jetFilter {
  enum { kJetChLowPt = 0,
         kJetChHighPt = 1,
         kTrackLowPt = 2,
         kTrackHighPt = 3,
         kTriggerObjects = 4
  };

  enum { kBinAllEvents = 0,
         kBinJetChLowPt = 1,
         kBinJetChHighPt = 2,
         kBinTrackLowPt = 3,
         kBinTrackHighPt = 4,
         kBins = 5
  };

  Produces<aod::JetFilters> tags;

  Configurable<std::string> evSel{"evSel", "sel8", "choose event selection"};
  Configurable<float> cfgJetR{"cfgJetR", 0.6,
                              "trigger jet resolution parameter"}; // jet cone radius

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "min track pT in filter"};
  Configurable<float> trackLowPtTriggerThreshold{"trackLowPtThreshold", 10.0, "low pT track trigger threshold"};
  Configurable<float> trackHighPtTriggerThreshold{"trackHighPtThreshold", 25.0, "high pT track trigger threshold"};
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

  Filter trackFilter = (nabs(aod::jtrack::eta) < static_cast<float>(cfgEtaTPC)) && (aod::jtrack::pt > trackPtMin);
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;

  void init(o2::framework::InitContext&)
  {
    triggerJetR = TMath::Nint(cfgJetR * 100.0f);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F,
                {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    spectra.add("ptphiJetChSelected_lowptjettrigger", "#it{p}_{T} of selected low pT charged jet trigger vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptphiJetChSelected_highptjettrigger", "#it{p}_{T} of selected high pT charged jet trigger vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptphiTrackSelected_lowtrackpttrigger", "#it{p}_{T} of selected low pT tracks vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptphiTrackSelected_hightrackpttrigger", "#it{p}_{T} of selected high pT tracks vs #varphi", HistType::kTH2F,
                {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi(), "#varphi"}});

    spectra.add("ptetaJetChSelected_lowptjettrigger", "#it{p}_{T} of selected low pT charged jet trigger vs #eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});

    spectra.add("ptetaJetChSelected_highptjettrigger", "#it{p}_{T} of selected high pT charged jet trigger vs #eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});

    spectra.add("ptetaTrackSelected_lowtrackpttrigger", "#it{p}_{T} of selected low pT tracks vs #eta", HistType::kTH2F,
                {{150, 0., +150., "track #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0, "#eta"}});

    spectra.add("ptetaTrackSelected_hightrackpttrigger", "#it{p}_{T} of selected high pT tracks vs #eta", HistType::kTH2F,
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
    spectra.add("hLeadingTrackPt", "Leading track pT in |#eta| < 0.9;",
                {HistType::kTH1F, {axisTrackPt}});

    spectra.add("hEtaVsPtTracksInclusive", "#eta of tracks |#eta| < 0.9;",
                {HistType::kTH2F, {axisTrackPt, axisEta}});

    spectra.add("hPhiVsPtTracksInclusive", "#varphi of tracks |#eta| < 0.9;",
                {HistType::kTH2F, {axisTrackPt, axisPhi}});

    spectra.add("hLeadingAKTJetR06Pt", "#it{p}_{T} of AKT R=0.6 charged jets in |#eta| < 0.9;",
                {HistType::kTH1F, {axisJetPt}});

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

    hProcessedEvents->GetXaxis()->SetBinLabel(kBinAllEvents + 1, "Processed events");
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinJetChLowPt + 1, o2::aod::filtering::JetChLowPt::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinJetChHighPt + 1, o2::aod::filtering::JetChHighPt::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinTrackLowPt + 1, o2::aod::filtering::TrackLowPt::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(kBinTrackHighPt + 1, o2::aod::filtering::TrackHighPt::columnLabel());
  }

  template <bool withRho, typename T, typename U, typename D>
  void doTriggering(T const& collision, U const& jets, D const& tracks)
  {

    // collision process loop
    bool keepEvent[kTriggerObjects]{false, false, false, false};
    if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>("NoTimeFrameBorder")))) {
      tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt], keepEvent[kTrackLowPt], keepEvent[kTrackHighPt]);
      return;
    }

    spectra.fill(HIST("fCollZpos"), collision.posZ());
    hProcessedEvents->Fill(static_cast<float>(kBinAllEvents) + 0.1f); // all minimum bias events

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt], keepEvent[kTrackLowPt], keepEvent[kTrackHighPt]);
      return;
    }

    // FILL SPECTRA OF INCLUSIVE JETS IN FIDUCIAL VOLUME
    if (TMath::Abs(collision.posZ()) < cfgZvtx) {
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
      spectra.fill(HIST("hLeadingAKTJetR06Pt"), jet.pt());

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

    float leadingTrackPt = -1.;
    float leadingTrackPhi = -100.;
    float leadingTrackEta = -100.;
    for (const auto& track : tracks) { // search for the leading track
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      spectra.fill(HIST("hPhiVsPtTracksInclusive"), track.pt(), track.phi());
      spectra.fill(HIST("hEtaVsPtTracksInclusive"), track.pt(), track.eta());
      if (track.pt() > leadingTrackPt) {
        leadingTrackPt = track.pt();
        leadingTrackPhi = track.phi();
        leadingTrackEta = track.eta();
      }
    }
    spectra.fill(HIST("hLeadingTrackPt"), leadingTrackPt);

    if (leadingTrackPt > trackLowPtTriggerThreshold) {
      keepEvent[kTrackLowPt] = true;
      spectra.fill(HIST("ptphiTrackSelected_lowtrackpttrigger"), leadingTrackPt, leadingTrackPhi);
      spectra.fill(HIST("ptetaTrackSelected_lowtrackpttrigger"), leadingTrackPt, leadingTrackEta);
    }
    if (leadingTrackPt > trackHighPtTriggerThreshold) {
      keepEvent[kTrackHighPt] = true;
      spectra.fill(HIST("ptphiTrackSelected_hightrackpttrigger"), leadingTrackPt, leadingTrackPhi);
      spectra.fill(HIST("ptetaTrackSelected_hightrackpttrigger"), leadingTrackPt, leadingTrackEta);
    }

    if (keepEvent[kJetChLowPt])
      hProcessedEvents->Fill(static_cast<float>(kBinJetChLowPt) + 0.1f);
    if (keepEvent[kJetChHighPt])
      hProcessedEvents->Fill(static_cast<float>(kBinJetChHighPt) + 0.1f);
    if (keepEvent[kTrackLowPt])
      hProcessedEvents->Fill(static_cast<float>(kBinTrackLowPt) + 0.1f);
    if (keepEvent[kTrackHighPt])
      hProcessedEvents->Fill(static_cast<float>(kBinTrackHighPt) + 0.1f);

    tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt], keepEvent[kTrackLowPt], keepEvent[kTrackHighPt]);
  }

  void processWithoutRho(aod::JetCollision const& collision, o2::aod::ChargedJets const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    doTriggering<false>(collision, jets, tracks);
  }
  PROCESS_SWITCH(jetFilter, processWithoutRho, "Do charged jet triggering without background estimation for filling histograms", true);

  void processWithRho(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision, o2::aod::ChargedJets const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    doTriggering<true>(collision, jets, tracks);
  }
  PROCESS_SWITCH(jetFilter, processWithRho, "Do charged jet triggering with background estimation for filling histograms", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<jetFilter>(cfg)};
}
