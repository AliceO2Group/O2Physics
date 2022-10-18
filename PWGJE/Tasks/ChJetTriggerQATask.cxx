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
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/HistogramRegistry.h"

#include <TMath.h>
#include <cmath>
#include <string>

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

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0,
                                   "Accepted z-vertex range"};
  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9,
                                   "Full eta range"}; // eta range of TPC
  Configurable<float> cfgJetR{"cfgJetR", 0.4,
                              "jet resolution parameter"}; // jet cone radius
  Configurable<float> cfgJetPtMin{
    "cfgJetPtMin", 0.05,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT
  Configurable<int> bTriggerDecision{
    "bTriggerDecision", 0,
    "Charged Jet Trigger Decision Selection"}; // 0=MB Event, 1=Event selected
                                               // by EPN
  double fiducialVolume;                       // 0.9 - jetR

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  HistogramRegistry spectra{
    "spectra",
    {
      {"vertexZ", "z vertex", {HistType::kTH1F, {{400, -20., +20.}}}}, //
      {"ptphiTrackInclGood",
       "pT vs phi inclusive good tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaTrackInclGood",
       "pT vs eta inclusive good tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"ptphiTrackInclRejected",
       "pT vs phi inclusive rejected tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaTrackInclRejected",
       "pT vs eta inclusive rejected tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"ptJetChPtInclFidVol",
       "inclusive charged jet pT in fiducial volume",
       {HistType::kTH1F, {{200, 0., +200.}}}}, //
      {"ptphiJetChPtInclFidVol",
       "inclusive charged jet pT vs phi in fiducial volume",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaJetChPtInclFidVol",
       "inclusive charged jet pT vs eta in fiducial volume",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"ptetaJetChPtInclFullVol",
       "inclusive charged jet pT vs eta in full TPC volume",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"fLeadJetChPtVsLeadingTrack",
       "inclusive charged jet pT in TPC volume",
       {HistType::kTH2F, {{200, 0., +200.}, {200, 0., +200.}}}}, //
      {"ptetaLeadingTrack",
       "pT vs eta leading tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"ptphiLeadingTrack",
       "pT vs phi leading tracks",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaLeadingJet",
       "pT vs eta leading jet",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"ptphiLeadingJet",
       "pT vs phi leading jet",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"fLeadJetEtaVsLeadingTrackEtaPathological",
       "leading jet eta versus track eta for cases pT jet is less than pT "
       "track in TPC volume",
       {HistType::kTH2F, {{40, -1., 1.}, {40, -1., 1.}}}}, //
      {"fLeadJetPhiVsLeadingTrackPhiPathological",
       "leading jet phi versus track phi for cases pT jet is less than pT "
       "track in TPC volume",
       {HistType::kTH2F,
        {{60, 0, TMath::TwoPi()}, {60, 0, TMath::TwoPi()}}}}, //
      {"trackTPCCrossedRowsAll",
       "trackTPCCrossedRowsAll",
       {HistType::kTH1F, {{200, 0., +200.}}}}, //
      {"trackTPCCrossedRowsGood",
       "trackTPCCrossedRowsGood",
       {HistType::kTH1F, {{200, 0., +200.}}}} //
    }};

  // TrackSelection globalTracks;
  void init(o2::framework::InitContext&)
  {

    // globalTracks = getGlobalTrackSelection();
    // globalTracks.SetEtaRange(-cfgTPCVolume, cfgTPCVolume);

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = cfgJetR;

    fiducialVolume = cfgTPCVolume - cfgJetR;
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,
                                    aod::TracksDCA, aod::TrackSelection>;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVertexCut);

  void
    process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                    aod::JetFilters>>::iterator const& collision,
            TrackCandidates const& tracks)
  {

    if (collision.hasJetChHighPt() >= bTriggerDecision) {
      jetConstituents.clear();
      jetReclustered.clear();

      double leadingJetPt = -1.0;
      double leadingJetEta = -2.0;
      double leadingJetPhi = -1.0;
      double leadingTrackPt = -1.0;
      double leadingTrackEta = -2.0;
      double leadingTrackPhi = -1.0;

      spectra.fill(HIST("vertexZ"),
                   collision.posZ()); // Inclusive Track Cross TPC Rows

      for (auto& trk : tracks) {
        if (fabs(trk.eta()) < cfgTPCVolume) {
          spectra.fill(
            HIST("trackTPCCrossedRowsAll"),
            trk.tpcNClsCrossedRows()); // Inclusive Track Cross TPC Rows

          if (!trk.isQualityTrack()) {

            spectra.fill(
              HIST("ptphiTrackInclRejected"), trk.pt(),
              trk.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
            spectra.fill(
              HIST("ptetaTrackInclRejected"), trk.pt(),
              trk.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

            continue; // skip bad quality tracks
          }

          spectra.fill(
            HIST("ptphiTrackInclGood"), trk.pt(),
            trk.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
          spectra.fill(
            HIST("ptetaTrackInclGood"), trk.pt(),
            trk.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

          spectra.fill(
            HIST("trackTPCCrossedRowsGood"),
            trk.tpcNClsCrossedRows()); // Inclusive Track Cross TPC Rows

          if (trk.pt() > cfgJetPtMin) { // jet constituents
            fillConstituents(
              trk,
              jetConstituents); // ./PWGJE/Core/JetFinder.h
                                // recombination scheme is assumed
                                // to be Escheme with pion mass
          }

          if (trk.pt() >
              leadingTrackPt) { // Find leading track pT in full TPC volume
            leadingTrackPt = trk.pt();
            leadingTrackEta = trk.eta();
            leadingTrackPhi = trk.phi();
          }
        }
      }

      if (leadingTrackPt > -1.) {
        spectra.fill(HIST("ptphiLeadingTrack"), leadingTrackPt,
                     leadingTrackPhi);
        spectra.fill(HIST("ptetaLeadingTrack"), leadingTrackPt,
                     leadingTrackEta);
      }

      // Reconstruct jet from tracks
      fastjet::ClusterSequenceArea clusterSeq(
        jetReclusterer.findJets(jetConstituents, jetReclustered));
      jetReclustered = sorted_by_pt(jetReclustered);

      // Find leading jet pT in full TPC volume
      for (auto& jet : jetReclustered) {
        if (fabs(jet.eta()) < cfgTPCVolume) {
          if (jet.pt() > leadingJetPt) {
            leadingJetPt = jet.pt();
            leadingJetEta = jet.eta();
            leadingJetPhi = jet.phi();
          }
        }
      }

      if (leadingJetPt > -1.) {
        spectra.fill(HIST("ptphiLeadingJet"), leadingJetPt, leadingJetPhi);
        spectra.fill(HIST("ptetaLeadingJet"), leadingJetPt, leadingJetEta);
      }

      spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt,
                   leadingJetPt); // leading jet pT versus leading track pT

      if ((leadingTrackPt - leadingJetPt) > 1e-5 && leadingJetPt > -1. &&
          leadingTrackPt > -1.) { // pathological case
        spectra.fill(HIST("fLeadJetEtaVsLeadingTrackEtaPathological"),
                     leadingTrackEta, leadingJetEta);
        spectra.fill(HIST("fLeadJetPhiVsLeadingTrackPhiPathological"),
                     leadingTrackPhi, leadingJetPhi);
      }
      //--------------------------------------------------------------

      // Inclusive Jet pT spectrum in Fiducial volume
      for (auto& jet : jetReclustered) {

        if (fabs(jet.eta()) < fiducialVolume) {
          spectra.fill(HIST("ptJetChPtInclFidVol"), jet.pt());
          spectra.fill(HIST("ptphiJetChPtInclFidVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChPtInclFidVol"), jet.pt(), jet.eta());
        }
        spectra.fill(HIST("ptetaJetChPtInclFullVol"), jet.pt(), jet.eta());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<ChJetTriggerQATask>(
    cfgc, TaskName{"jet-charged-trigger-qa"})};
}
