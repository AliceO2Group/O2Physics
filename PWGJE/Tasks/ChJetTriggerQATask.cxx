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
/// \author Filip Krizek <Filip.Krizek@cern.ch>
#include <cmath>
#include <string>
#include <vector>
#include <TMath.h>
#include <TVector2.h>
#include <TLorentzVector.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "EventFiltering/filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/HistogramRegistry.h"

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
// It would be good to run it for several jet radii  e.g. 0.2, 0.4, 0.6

struct ChJetTriggerQATask {

  Configurable<std::string> evSel{"evSel", "sel8", "choose event selection"};
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0,
                                   "Accepted z-vertex range"};
  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9,
                                   "Full eta range"}; // eta range of TPC
  Configurable<float> cfgJetR{"cfgJetR", 0.4,
                              "jet resolution parameter"}; // jet cone radius
  Configurable<float> cfgJetPtMin{
    "cfgJetPtMin", 0.15,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT

  Configurable<float> cfgTrackPhiMinCut{"cfgTrackPhiMinCut", -999, "track min phi cut"};
  Configurable<float> cfgTrackPhiMaxCut{"cfgTrackPhiMaxCut", 999, "track max phi cut"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<bool> bLowPtTrigger{"bLowPtTrigger", false, "charged jet low pT trigger selection"};
  Configurable<bool> bHighPtTrigger{"bHighPtTrigger", false, "charged jet high pT trigger selection"};
  Configurable<bool> bTrackPtTrigger{"bTrackPtTrigger", false, "track pT trigger selection"};

  Configurable<bool> bAddSupplementHistosToOutput{"bAddAdditionalHistosToOutput", false, "add supplementary histos to the output"};

  float fiducialVolume; // 0.9 - jetR

  HistogramRegistry spectra;

  int eventSelection = -1;
  int trackSelection = -1;
  void init(o2::framework::InitContext&)
  {
    fiducialVolume = static_cast<float>(cfgTPCVolume) - static_cast<float>(cfgJetR);
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(evSel));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // Basic histos
    spectra.add("vertexZ", "z vertex", {HistType::kTH1F, {{400, -20., +20.}}});
    spectra.add("ptphiTrackInclGood", "pT vs phi inclusive good tracks", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});
    spectra.add("ptetaTrackInclGood", "pT vs eta inclusive good tracks", {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}});
    spectra.add("ptLeadingTrack", "pT leading track", {HistType::kTH1F, {{100, 0., +100.}}});
    spectra.add("ptJetChInclFidVol", "inclusive charged jet pT in fiducial volume", {HistType::kTH1F, {{200, 0., +200.}}});
    spectra.add("ptphiJetChInclFidVol", "inclusive charged jet pT vs phi in fiducial volume", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});
    spectra.add("ptphiJetChInclFullVol", "inclusive charged jet pT vs phi in full TPC volume", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});
    spectra.add("ptetaJetChInclFidVol", "inclusive charged jet pT vs eta in fiducial volume", {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}});
    spectra.add("ptetaJetChInclFullVol", "inclusive charged jet pT vs eta in full TPC volume", {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}});
    spectra.add("ptetaLeadingJetFullVol", "pT vs eta leading jet", {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}});
    spectra.add("ptphiLeadingJetFullVol", "pT vs phi leading jet", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});

    // Supplementary plots
    if (bAddSupplementHistosToOutput) {
      spectra.add("ptJetChInclFullVol", "inclusive charged jet pT in full volume", {HistType::kTH1F, {{200, 0., +200.}}});
      spectra.add("phietaTrackAllInclGood", "phi vs eta all inclusive good tracks", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("phietaTrackHighPtInclGood", "phi vs eta inclusive good tracks with pT > 10 GeV", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("phietaJetChInclFidVol", "inclusive charged jet phi vs eta in fiducial volume", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("phietaJetChInclFullVol", "inclusive charged jet phi vs eta in full TPC volume", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("phietaJetChInclHighPtFidVol", "inclusive charged jet phi vs eta in fiducial volume", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("phietaJetChInclHighPtFullVol", "inclusive charged jet phi vs eta in full TPC volume", {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("ptetaLeadingTrack", "pT vs eta leading tracks", {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}});
      spectra.add("ptphiLeadingTrack", "pT vs phi leading tracks", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});
      spectra.add("jetAreaFullVol", "area of all jets in full TPC volume", {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}});
      spectra.add("jetAreaFidVol", "area of all jets in fiducial volume", {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}});
      spectra.add("fLeadJetChPtVsLeadingTrack", "inclusive charged jet pT in TPC volume", {HistType::kTH2F, {{100, 0., +100.}, {100, 0., +100.}}});
    }
  }

  // declare filters on collisions
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < static_cast<float>(cfgVertexCut));

  // declare filters on tracks
  Filter trackFilter = (nabs(aod::jtrack::eta) < static_cast<float>(cfgTPCVolume)) && (aod::jtrack::phi > static_cast<float>(cfgTrackPhiMinCut)) && (aod::jtrack::phi < static_cast<float>(cfgTrackPhiMaxCut)) && (aod::jtrack::pt > static_cast<float>(cfgJetPtMin));

  // declare filters on jets
  Filter jetRadiusSelection = (o2::aod::jet::r == nround(cfgJetR.node() * 100.0f));

  using filteredJets = o2::soa::Filtered<o2::aod::ChargedJets>;

  void
    process(soa::Filtered<soa::Join<JetCollisions,
                                    aod::JChTrigSels, aod::EvSels>>::iterator const& collision,
            soa::Filtered<JetTracks> const& tracks, o2::soa::Filtered<soa::Join<o2::aod::ChargedJets, aod::ChargedJetConstituents>> const& jets)
  {

    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    if ((bLowPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) || (bHighPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) || (bTrackPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::trackPt)) || ((!bLowPtTrigger) && (!bHighPtTrigger) && (!bTrackPtTrigger))) {
      // bLowPtTrigger=1  and bHighPtTrigger=0 --> fill histos with low trigger only
      // bLowPtTrigger=0  and bHighPtTrigger=1 --> fill histos with high trigger only
      // bLowPtTrigger=1  and bHighPtTrigger=1 --> fill histos with mixture of low and high trigger
      // bTrackPtTrigger=1 --> fill histos for high pt track trigger
      // bLowPtTrigger=0  and bHighPtTrigger=0 and bTrackPtTrigger=0 --> fill histos with minimum bias ie. ignore trigger decision

      float leadingJetPt = -1.0;
      float leadingJetEta = -2.0;
      float leadingJetPhi = -1.0;
      float leadingTrackPt = -1.0;
      float leadingTrackEta = -2.0;
      float leadingTrackPhi = -1.0;

      spectra.fill(HIST("vertexZ"),
                   collision.posZ()); // Inclusive Track Cross TPC Rows

      for (auto& trk : tracks) { // loop over filtered tracks in full TPC volume having pT > 100 MeV

        if (!jetderiveddatautilities::selectTrack(trk, trackSelection)) {
          continue;
        }

        spectra.fill(
          HIST("ptphiTrackInclGood"), trk.pt(),
          trk.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
        spectra.fill(
          HIST("ptetaTrackInclGood"), trk.pt(),
          trk.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

        if (bAddSupplementHistosToOutput) {
          spectra.fill(
            HIST("phietaTrackAllInclGood"), trk.eta(),
            trk.phi()); // Inclusive Track pT vs eta spectrum in TPC volume

          if (trk.pt() > 5.0) {
            spectra.fill(
              HIST("phietaTrackHighPtInclGood"), trk.eta(),
              trk.phi()); // Inclusive Track pT vs eta spectrum in TPC volume
          }
        }

        if (trk.pt() >
            leadingTrackPt) { // Find leading track pT in full TPC volume
          leadingTrackPt = trk.pt();
          leadingTrackEta = trk.eta();
          leadingTrackPhi = trk.phi();
        }
      }

      if (leadingTrackPt > -1.) {
        spectra.fill(HIST("ptLeadingTrack"), leadingTrackPt);
      }

      if (bAddSupplementHistosToOutput) {
        if (leadingTrackPt > -1.) {
          spectra.fill(HIST("ptphiLeadingTrack"), leadingTrackPt,
                       leadingTrackPhi);
          spectra.fill(HIST("ptetaLeadingTrack"), leadingTrackPt,
                       leadingTrackEta);
        }
      }

      // Find leading jet pT in full TPC volume
      for (auto& jet : jets) {
        if (fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {

          if (jet.pt() > leadingJetPt) {
            leadingJetPt = jet.pt();
            leadingJetEta = jet.eta();
            leadingJetPhi = jet.phi();
          }
        }
      }

      if (leadingJetPt > -1.) {
        spectra.fill(HIST("ptphiLeadingJetFullVol"), leadingJetPt, leadingJetPhi);
        spectra.fill(HIST("ptetaLeadingJetFullVol"), leadingJetPt, leadingJetEta);
      }

      if (bAddSupplementHistosToOutput) {
        if (leadingJetPt > -1. && leadingTrackPt > -1.) {
          spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt,
                       leadingJetPt); // leading jet pT versus leading track pT
        }
      }

      // Inclusive Jet pT spectrum in Fiducial volume
      for (auto& jet : jets) {
        if (fabs(jet.eta()) < fiducialVolume) {
          spectra.fill(HIST("ptJetChInclFidVol"), jet.pt());
          spectra.fill(HIST("ptphiJetChInclFidVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFidVol"), jet.pt(), jet.eta());

          if (bAddSupplementHistosToOutput) {
            spectra.fill(HIST("phietaJetChInclFidVol"), jet.eta(), jet.phi());
            if (jet.pt() > 10.0) {
              spectra.fill(HIST("phietaJetChInclHighPtFidVol"), jet.eta(), jet.phi());
            }
            spectra.fill(HIST("jetAreaFidVol"), jet.pt(), jet.area());
          }
        }

        if (fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {
          spectra.fill(HIST("ptphiJetChInclFullVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFullVol"), jet.pt(), jet.eta());

          if (bAddSupplementHistosToOutput) {
            spectra.fill(HIST("ptJetChInclFullVol"), jet.pt());

            spectra.fill(HIST("phietaJetChInclFullVol"), jet.eta(), jet.phi());
            if (jet.pt() > 10.0) {
              spectra.fill(HIST("phietaJetChInclHighPtFullVol"), jet.eta(), jet.phi());
            }
            spectra.fill(HIST("jetAreaFullVol"), jet.pt(), jet.area());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<ChJetTriggerQATask>(
    cfgc, TaskName{"jet-charged-trigger-qa"})};
}
