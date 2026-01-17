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

/// \author Filip Krizek <filip.krizek@cern.ch>
/// \author Kotliarov Artem <artem.kotliarov@cern.ch>
/// \file jetTriggerChargedQa.cxx
/// \brief QA of trigger performance for charged jets

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FilteredColl = soa::Filtered<soa::Join<aod::JetCollisions, aod::JChTrigSels>>::iterator;
using FilteredJTracks = soa::Filtered<soa::Join<aod::JTracks, aod::JTrackPIs, aod::JTrackExtras>>;
using FilteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using JoinedTracks = soa::Join<aod::Tracks, aod::TracksExtra>;

float dcaXYPtCut(float tracPt)
{
  return 0.0105f + 0.0350f / std::pow(tracPt, 1.1f);
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

struct JetTriggerChargedQa {

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
  Configurable<bool> bAddSupplementHistosToOutput{"bAddSupplementHistosToOutput", false, "add supplementary histos to the output"};
  Configurable<bool> bStudyPhiTrack{"bStudyPhiTrack", false, "add histos for detailed study of track phi distribution"};

  Configurable<float> phiAngleRestriction{"phiAngleRestriction", 0.3, "angle to restrict track phi for plotting tpc momentum"};
  Configurable<float> dcaXYMultFact{"dcaXYMultFact", 3., "mult factor to relax pT dependent dcaXY cut for quality tracks"};
  Configurable<float> dcaZCut{"dcaZCut", 3., "cut on dcaZ for quality tracks"};

  float twoPi = constants::math::TwoPI;
  ConfigurableAxis dcaXYBinning{"dcaXYBinning", {100, -5., 5.}, ""};
  ConfigurableAxis dcaZBinning{"dcaZBinning", {100, -3., 3.}, ""};

  ConfigurableAxis xPhiAxis{"xPhiAxis", {40, 0., twoPi}, ""};
  ConfigurableAxis yQ1pTAxis{"yQ1pTAxis", {200, -0.5, 0.5}, ""};

  float fiducialVolume = 0.0; // 0.9 - jetR

  HistogramRegistry spectra;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(InitContext&)
  {
    fiducialVolume = static_cast<float>(cfgTPCVolume) - static_cast<float>(cfgJetR);
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // Basic histos
    spectra.add("vertexZ", "z vertex", kTH1F, {{60, -12., 12.}});
    spectra.add("ptphiTrackInclGood", "pT vs phi inclusive good tracks", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
    spectra.add("ptetaTrackInclGood", "pT vs eta inclusive good tracks", kTH2F, {{100, 0., 100.}, {40, -1., 1.}});
    spectra.add("ptLeadingTrack", "pT leading track", kTH1F, {{100, 0., 100.}});
    spectra.add("ptJetChInclFidVol", "inclusive charged jet pT in fiducial volume", kTH1F, {{200, 0., 200.}});
    spectra.add("ptphiJetChInclFidVol", "inclusive charged jet pT vs phi in fiducial volume", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
    spectra.add("ptphiJetChInclFullVol", "inclusive charged jet pT vs phi in full TPC volume", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
    spectra.add("ptetaJetChInclFidVol", "inclusive charged jet pT vs eta in fiducial volume", kTH2F, {{100, 0., 100.}, {40, -1., 1.}});
    spectra.add("ptetaJetChInclFullVol", "inclusive charged jet pT vs eta in full TPC volume", kTH2F, {{100, 0., 100.}, {40, -1., 1.}});
    spectra.add("ptetaLeadingJetFullVol", "pT vs eta leading jet", kTH2F, {{100, 0., 100.}, {40, -1., 1.}});
    spectra.add("ptphiLeadingJetFullVol", "pT vs phi leading jet", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});

    // Supplementary plots
    if (bAddSupplementHistosToOutput) {
      spectra.add("ptJetChInclFullVol", "inclusive charged jet pT in full volume", kTH1F, {{200, 0., 200.}});
      spectra.add("phietaTrackAllInclGood", "phi vs eta all inclusive good tracks", kTH2F, {{80, -1., 1.}, {40, 0, twoPi}});
      spectra.add("phietaTrackHighPtInclGood", "phi vs eta inclusive good tracks with pT > 10 GeV", kTH2F, {{40, -1., 1.}, {40, 0, twoPi}});
      spectra.add("phietaJetChInclFidVol", "inclusive charged jet phi vs eta in fiducial volume", kTH2F, {{40, -1., 1.}, {40, 0, twoPi}});
      spectra.add("phietaJetChInclFullVol", "inclusive charged jet phi vs eta in full TPC volume", kTH2F, {{40, -1., 1.}, {40, 0, twoPi}});
      spectra.add("phietaJetChInclHighPtFidVol", "inclusive charged jet phi vs eta in fiducial volume", kTH2F, {{40, -1., 1.}, {40, 0, twoPi}});
      spectra.add("phietaJetChInclHighPtFullVol", "inclusive charged jet phi vs eta in full TPC volume", kTH2F, {{40, -1., 1.}, {40, 0, twoPi}});
      spectra.add("ptetaLeadingTrack", "pT vs eta leading tracks", kTH2F, {{100, 0., 100.}, {40, -1., 1.}});
      spectra.add("ptphiLeadingTrack", "pT vs phi leading tracks", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
      spectra.add("jetAreaFullVol", "area of all jets in full TPC volume", kTH2F, {{100, 0., 100.}, {50, 0., 2.}});
      spectra.add("jetAreaFidVol", "area of all jets in fiducial volume", kTH2F, {{100, 0., 100.}, {50, 0., 2.}});
      spectra.add("fLeadJetChPtVsLeadingTrack", "inclusive charged jet pT in TPC volume", kTH2F, {{100, 0., 100.}, {100, 0., 100.}});
    }

    // Study of non-uniformity of phi distribution of tracks
    if (bStudyPhiTrack) {
      spectra.add("globalP_tpcglobalPDiff_withoutcuts", "difference of global and TPC inner momentum vs global momentum without any selection applied", kTH2F, {{100, 0., 100.}, {200, -100., 100.}});
      spectra.add("globalP_tpcglobalPDiff", "difference of global and TPC inner momentum vs global momentum with selection applied", kTH2F, {{100, 0., 100.}, {200, -100., 100.}});
      spectra.add("global1overP_tpcglobalPDiff", "difference of global and TPC inner momentum vs global momentum with selection applied", kTH2F, {{100, 0., 100.}, {125, -8., 8.}});

      spectra.add("globalP_tpcglobalPDiff_withoutcuts_phirestrict", "difference of global and TPC inner momentum vs global momentum without any selection applied in a restricted phi", kTH2F, {{100, 0., 100.}, {200, -100., 100.}});
      spectra.add("globalP_tpcglobalPDiff_phirestrict", "difference of global and TPC inner momentum vs global momentum with selection applied restricted phi", kTH2F, {{100, 0., 100.}, {200, -100., 100.}});
      spectra.add("global1overP_tpcglobalPDiff_phirestrict", "difference of 1/p global and TPC inner momentum vs global momentum with selection applied restricted phi", kTH2F, {{100, 0., 100.}, {500, -8., 8.}});

      spectra.add("DCAxy_track_Phi_pT", "track DCAxy vs phi & pT of tracks w. nITSClusters #geq 4", kTH3F, {dcaXYBinning, {40, 0., twoPi}, {100, 0., 100.}});
      spectra.add("DCAz_track_Phi_pT", "track DCAz vs phi & pT of tracks w. nITSClusters #geq 4", kTH3F, {dcaZBinning, {40, 0., twoPi}, {100, 0., 100.}});
      spectra.add("nITSClusters_TrackPt", "Number of ITS hits vs phi & pT of tracks", kTH3F, {{7, 1., 8.}, {40, 0., twoPi}, {100, 0., 100.}});
      spectra.add("ptphiQualityTracks", "pT vs phi of quality tracks", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
      spectra.add("ptphiAllTracks", "pT vs phi of all tracks", kTH2F, {{100, 0., 100.}, {40, 0, twoPi}});
      spectra.add("phi_Q1pT", "Track phi vs. q/pT", kTH2F, {xPhiAxis, yQ1pTAxis});
    }
  }

  // declare filters on collisions
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < static_cast<float>(cfgVertexCut));

  // declare filters on tracks
  Filter trackFilter = (nabs(aod::jtrack::eta) < static_cast<float>(cfgTPCVolume)) && (aod::jtrack::phi > static_cast<float>(cfgTrackPhiMinCut)) && (aod::jtrack::phi < static_cast<float>(cfgTrackPhiMaxCut)) && (aod::jtrack::pt > static_cast<float>(cfgJetPtMin));

  // declare filters on jets
  Filter jetRadiusSelection = (aod::jet::r == nround(cfgJetR.node() * 100.0f));

  void process(FilteredColl const& collision, FilteredJTracks const& tracks, FilteredJets const& jets, JoinedTracks const&)
  {

    if (!jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>("NoTimeFrameBorder")))) {
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

      spectra.fill(HIST("vertexZ"), collision.posZ()); // Inclusive Track Cross TPC Rows

      // loop over filtered tracks in full TPC volume having pT > 100 MeV
      for (auto const& track : tracks) {

        auto const& originalTrack = track.track_as<JoinedTracks>();

        if (bStudyPhiTrack) {

          spectra.fill(HIST("globalP_tpcglobalPDiff_withoutcuts"), track.p(), track.p() - originalTrack.tpcInnerParam());
          spectra.fill(HIST("ptphiAllTracks"), track.pt(), track.phi());

          if (std::fabs(track.phi() - constants::math::PI) < phiAngleRestriction) {
            spectra.fill(HIST("globalP_tpcglobalPDiff_withoutcuts_phirestrict"), track.p(), track.p() - originalTrack.tpcInnerParam());
          }
        }

        if (!jetderiveddatautilities::selectTrack(track, trackSelection))
          continue;

        spectra.fill(HIST("ptphiTrackInclGood"), track.pt(), track.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
        spectra.fill(HIST("ptetaTrackInclGood"), track.pt(), track.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

        if (bAddSupplementHistosToOutput) {
          spectra.fill(HIST("phietaTrackAllInclGood"), track.eta(), track.phi()); // Inclusive Track pT vs eta spectrum in TPC volume

          float trackPtCut = 5.0;
          if (track.pt() > trackPtCut) {
            spectra.fill(HIST("phietaTrackHighPtInclGood"), track.eta(), track.phi()); // Inclusive Track pT vs eta spectrum in TPC volume
          }
        }

        if (track.pt() > leadingTrackPt) { // Find leading track pT in full TPC volume
          leadingTrackPt = track.pt();
          leadingTrackEta = track.eta();
          leadingTrackPhi = track.phi();
        }

        if (bStudyPhiTrack) {
          spectra.fill(HIST("phi_Q1pT"), originalTrack.phi(), originalTrack.sign() / originalTrack.pt());
          spectra.fill(HIST("ptphiQualityTracks"), track.pt(), track.phi());

          bool bDcaCondition = (std::fabs(track.dcaZ()) < dcaZCut) && (std::fabs(track.dcaXY()) < dcaXYMultFact * dcaXYPtCut(track.pt()));

          int nITSClusters = 4;
          if (originalTrack.itsNCls() >= nITSClusters && bDcaCondition) { // correspond to number of track hits in ITS layers
            spectra.fill(HIST("DCAxy_track_Phi_pT"), track.dcaXY(), track.phi(), track.pt());
            spectra.fill(HIST("DCAz_track_Phi_pT"), track.dcaZ(), track.phi(), track.pt());
          }

          spectra.fill(HIST("nITSClusters_TrackPt"), originalTrack.itsNCls(), track.phi(), track.pt());

          spectra.fill(HIST("globalP_tpcglobalPDiff"), track.p(), track.p() - originalTrack.tpcInnerParam());
          if (track.p() > 0 && originalTrack.tpcInnerParam() > 0) {
            spectra.fill(HIST("global1overP_tpcglobalPDiff"), track.p(), 1. / track.p() - 1. / originalTrack.tpcInnerParam());
          }

          if (std::fabs(track.phi() - constants::math::PI) < phiAngleRestriction) {
            spectra.fill(HIST("globalP_tpcglobalPDiff_phirestrict"), track.p(), track.p() - originalTrack.tpcInnerParam());

            if (track.p() > 0 && originalTrack.tpcInnerParam() > 0) {
              spectra.fill(HIST("global1overP_tpcglobalPDiff_phirestrict"), track.p(), 1. / track.p() - 1. / originalTrack.tpcInnerParam());
            }
          }
        }
      }

      if (leadingTrackPt > -1.) {
        spectra.fill(HIST("ptLeadingTrack"), leadingTrackPt);
      }

      if (bAddSupplementHistosToOutput) {
        if (leadingTrackPt > -1.) {
          spectra.fill(HIST("ptphiLeadingTrack"), leadingTrackPt, leadingTrackPhi);
          spectra.fill(HIST("ptetaLeadingTrack"), leadingTrackPt, leadingTrackEta);
        }
      }

      // Find leading jet pT in full TPC volume
      for (const auto& jet : jets) {
        if (std::fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {

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
          spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt, leadingJetPt); // leading jet pT versus leading track pT
        }
      }

      // Inclusive Jet pT spectrum in Fiducial volume
      for (const auto& jet : jets) {
        if (std::fabs(jet.eta()) < fiducialVolume) {
          spectra.fill(HIST("ptJetChInclFidVol"), jet.pt());
          spectra.fill(HIST("ptphiJetChInclFidVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFidVol"), jet.pt(), jet.eta());

          if (bAddSupplementHistosToOutput) {
            spectra.fill(HIST("phietaJetChInclFidVol"), jet.eta(), jet.phi());

            float jetPtCut = 10.0;
            if (jet.pt() > jetPtCut) {
              spectra.fill(HIST("phietaJetChInclHighPtFidVol"), jet.eta(), jet.phi());
            }
            spectra.fill(HIST("jetAreaFidVol"), jet.pt(), jet.area());
          }
        }

        if (std::fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {
          spectra.fill(HIST("ptphiJetChInclFullVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFullVol"), jet.pt(), jet.eta());

          if (bAddSupplementHistosToOutput) {
            spectra.fill(HIST("ptJetChInclFullVol"), jet.pt());

            spectra.fill(HIST("phietaJetChInclFullVol"), jet.eta(), jet.phi());

            float jetPtCut = 10.0;
            if (jet.pt() > jetPtCut) {
              spectra.fill(HIST("phietaJetChInclHighPtFullVol"), jet.eta(), jet.phi());
            }
            spectra.fill(HIST("jetAreaFullVol"), jet.pt(), jet.area());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetTriggerChargedQa>(cfgc)}; }
