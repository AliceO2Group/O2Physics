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

using filteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::JChTrigSels, aod::EvSels>>::iterator;
using filteredJTracks = soa::Filtered<soa::Join<aod::JTracks, aod::JTrackPIs, aod::JTrackExtras>>;
using filteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using joinedTracks = soa::Join<aod::Tracks, aod::TracksExtra>;

float DcaXYPtCut(float tracPt)
{
  return 0.0105f + 0.0350f / pow(tracPt, 1.1f);
}

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
  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0, "Accepted z-vertex range"};
  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9, "Full eta range"};                // eta range of TPC
  Configurable<float> cfgJetR{"cfgJetR", 0.4, "jet resolution parameter"};                // jet cone radius
  Configurable<float> cfgJetPtMin{"cfgJetPtMin", 0.15, "minimum jet pT constituent cut"}; // minimum jet constituent pT

  Configurable<float> cfgTrackPhiMinCut{"cfgTrackPhiMinCut", -999, "track min phi cut"};
  Configurable<float> cfgTrackPhiMaxCut{"cfgTrackPhiMaxCut", 999, "track max phi cut"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<bool> bLowPtTrigger{"bLowPtTrigger", false, "charged jet low pT trigger selection"};
  Configurable<bool> bHighPtTrigger{"bHighPtTrigger", false, "charged jet high pT trigger selection"};
  Configurable<bool> bTrackLowPtTrigger{"bTrackLowPtTrigger", false, "track low pT trigger selection"};
  Configurable<bool> bTrackHighPtTrigger{"bTrackHighPtTrigger", false, "track high pT trigger selection"};
  Configurable<bool> bAddSupplementHistosToOutput{"bAddAdditionalHistosToOutput", false, "add supplementary histos to the output"};

  Configurable<float> phiAngleRestriction{"phiAngleRestriction", 0.3, "angle to restrict track phi for plotting tpc momentum"};
  Configurable<float> dcaXY_multFact{"dcaXY_multFact", 3., "mult factor to relax pT dependent dcaXY cut for quality tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 3., "cut on dcaZ for quality tracks"};

  ConfigurableAxis dcaXY_Binning{"dcaXY_Binning", {100, -5., 5.}, ""};
  ConfigurableAxis dcaZ_Binning{"dcaZ_Binning", {100, -3., 3.}, ""};

  ConfigurableAxis xPhiAxis{"xPhiAxis", {180, 0., TMath::TwoPi()}, ""};
  ConfigurableAxis yQ1pTAxis{"yQ1pTAxis", {200, -0.5, 0.5}, ""};

  float fiducialVolume; // 0.9 - jetR

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(InitContext&)
  {
    fiducialVolume = static_cast<float>(cfgTPCVolume) - static_cast<float>(cfgJetR);
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));
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

    spectra.add("globalP_tpcglobalPDiff_withoutcuts", "difference of global and TPC inner momentum vs global momentum without any selection applied", {HistType::kTH2F, {{100, 0., +100.}, {200, -100., +100.}}});
    spectra.add("globalP_tpcglobalPDiff", "difference of global and TPC inner momentum vs global momentum with selection applied", {HistType::kTH2F, {{100, 0., +100.}, {200, -100., +100.}}});
    spectra.add("global1overP_tpcglobalPDiff", "difference of global and TPC inner momentum vs global momentum with selection applied", {HistType::kTH2F, {{100, 0., +100.}, {500, -8., +8.}}});

    spectra.add("globalP_tpcglobalPDiff_withoutcuts_phirestrict", "difference of global and TPC inner momentum vs global momentum without any selection applied in a restricted phi", {HistType::kTH2F, {{100, 0., +100.}, {200, -100., +100.}}});
    spectra.add("globalP_tpcglobalPDiff_phirestrict", "difference of global and TPC inner momentum vs global momentum with selection applied restricted phi", {HistType::kTH2F, {{100, 0., +100.}, {200, -100., +100.}}});
    spectra.add("global1overP_tpcglobalPDiff_phirestrict", "difference of 1/p global and TPC inner momentum vs global momentum with selection applied restricted phi", {HistType::kTH2F, {{100, 0., +100.}, {500, -8., +8.}}});

    spectra.add("DCAxy_track_Phi_pT", "track DCAxy vs phi & pT of tracks w. nITSClusters #geq 4", kTH3F, {dcaXY_Binning, {60, 0., TMath::TwoPi()}, {100, 0., 100.}});
    spectra.add("DCAz_track_Phi_pT", "track DCAz vs phi & pT of tracks w. nITSClusters #geq 4", kTH3F, {dcaZ_Binning, {60, 0., TMath::TwoPi()}, {100, 0., 100.}});
    spectra.add("nITSClusters_TrackPt", "Number of ITS hits vs phi & pT of tracks", kTH3F, {{7, 1., 8.}, {60, 0., TMath::TwoPi()}, {100, 0., 100.}});
    spectra.add("ptphiQualityTracks", "pT vs phi of quality tracks", {HistType::kTH2F, {{100, 0., 100.}, {60, 0, TMath::TwoPi()}}});
    spectra.add("ptphiAllTracks", "pT vs phi of all tracks", {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}});
    spectra.add("phi_Q1pT", "Track phi vs. q/pT", kTH2F, {xPhiAxis, yQ1pTAxis});

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
  Filter jetRadiusSelection = (aod::jet::r == nround(cfgJetR.node() * 100.0f));

  void process(filteredColl const& collision, filteredJTracks const& tracks, filteredJets const& jets, joinedTracks const&)
  {

    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    bool bLowPtJet = (bLowPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::jetChLowPt));
    bool bHighPtJet = (bHighPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::jetChHighPt));
    bool bLowPtTrack = (bTrackLowPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::trackLowPt));
    bool bHighPtTrack = (bTrackHighPtTrigger && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::trackHighPt));
    bool bMinimumBias = ((!bLowPtTrigger) && (!bHighPtTrigger) && (!bTrackLowPtTrigger) && (!bTrackHighPtTrigger));

    if (bLowPtJet || bHighPtJet || bLowPtTrack || bHighPtTrack || bMinimumBias) {
      // bLowPtTrigger=1  and bHighPtTrigger=0 --> fill histos with low trigger only
      // bLowPtTrigger=0  and bHighPtTrigger=1 --> fill histos with high trigger only
      // bLowPtTrigger=1  and bHighPtTrigger=1 --> fill histos with mixture of low and high trigger
      // bTrackLowPtTrigger=1 --> fill histos for low pt track trigger
      // bTrackHighPtTrigger=1 --> fill histos for high pt track trigger
      // bLowPtTrigger=0 and bHighPtTrigger=0 and bTrackLowPtTrigger=0 and bTrackHighPtTrigger=0 --> fill histos with minimum bias ie. ignore trigger decision

      float leadingJetPt = -1.0;
      float leadingJetEta = -2.0;
      float leadingJetPhi = -1.0;
      float leadingTrackPt = -1.0;
      float leadingTrackEta = -2.0;
      float leadingTrackPhi = -1.0;

      spectra.fill(HIST("vertexZ"),
                   collision.posZ()); // Inclusive Track Cross TPC Rows

      for (auto const& track : tracks) { // loop over filtered tracks in full TPC volume having pT > 100 MeV

        auto const& originalTrack = track.track_as<joinedTracks>();

        spectra.fill(HIST("globalP_tpcglobalPDiff_withoutcuts"), track.p(), track.p() - originalTrack.tpcInnerParam());

        if (TMath::Abs(track.phi() - TMath::Pi()) < phiAngleRestriction) {
          spectra.fill(HIST("globalP_tpcglobalPDiff_withoutcuts_phirestrict"), track.p(), track.p() - originalTrack.tpcInnerParam());
        }

        spectra.fill(HIST("ptphiAllTracks"), track.pt(), track.phi());

        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }

        spectra.fill(HIST("phi_Q1pT"), originalTrack.phi(), originalTrack.sign() / originalTrack.pt());
        spectra.fill(HIST("ptphiQualityTracks"), track.pt(), track.phi());

        bool bDcaCondition = (fabs(track.dcaZ()) < dcaZ_cut) && (fabs(track.dcaXY()) < dcaXY_multFact * DcaXYPtCut(track.pt()));
        if (originalTrack.itsNCls() >= 4 && bDcaCondition) { // correspond to number of track hits in ITS layers
          spectra.fill(HIST("DCAxy_track_Phi_pT"), track.dcaXY(), track.phi(), track.pt());
          spectra.fill(HIST("DCAz_track_Phi_pT"), track.dcaZ(), track.phi(), track.pt());
        }

        spectra.fill(HIST("nITSClusters_TrackPt"), originalTrack.itsNCls(), track.phi(), track.pt());

        spectra.fill(HIST("globalP_tpcglobalPDiff"), track.p(), track.p() - originalTrack.tpcInnerParam());
        if (track.p() > 0 && originalTrack.tpcInnerParam() > 0) {
          spectra.fill(HIST("global1overP_tpcglobalPDiff"), track.p(), 1. / track.p() - 1. / originalTrack.tpcInnerParam());
        }
        if (TMath::Abs(track.phi() - TMath::Pi()) < phiAngleRestriction) {
          spectra.fill(HIST("globalP_tpcglobalPDiff_phirestrict"), track.p(), track.p() - originalTrack.tpcInnerParam());

          if (track.p() > 0 && originalTrack.tpcInnerParam() > 0) {
            spectra.fill(HIST("global1overP_tpcglobalPDiff_phirestrict"), track.p(), 1. / track.p() - 1. / originalTrack.tpcInnerParam());
          }
        }

        spectra.fill(
          HIST("ptphiTrackInclGood"), track.pt(),
          track.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
        spectra.fill(
          HIST("ptetaTrackInclGood"), track.pt(),
          track.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

        if (bAddSupplementHistosToOutput) {
          spectra.fill(
            HIST("phietaTrackAllInclGood"), track.eta(),
            track.phi()); // Inclusive Track pT vs eta spectrum in TPC volume

          if (track.pt() > 5.0) {
            spectra.fill(
              HIST("phietaTrackHighPtInclGood"), track.eta(),
              track.phi()); // Inclusive Track pT vs eta spectrum in TPC volume
          }
        }

        if (track.pt() >
            leadingTrackPt) { // Find leading track pT in full TPC volume
          leadingTrackPt = track.pt();
          leadingTrackEta = track.eta();
          leadingTrackPhi = track.phi();
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
