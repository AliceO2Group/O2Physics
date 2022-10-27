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
#include <TVector2.h>
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
// It would be good to run it for several jet radii  e.g. 0.2, 0.4, 0.6

struct ChJetTriggerQATask {

  Configurable<float> cfgVertexCut{"cfgVertexCut", 10.0,
                                   "Accepted z-vertex range"};
  Configurable<float> cfgTPCVolume{"cfgTPCVolume", 0.9,
                                   "Full eta range"}; // eta range of TPC
  Configurable<float> cfgJetR{"cfgJetR", 0.4,
                              "jet resolution parameter"}; // jet cone radius
  Configurable<float> cfgJetPtMin{
    "cfgJetPtMin", 0.1,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT
  Configurable<int> bTriggerDecision{
    "bTriggerDecision", 0,
    "Charged Jet Trigger Decision Selection"}; // 0=MB Event, 1=Event selected
                                               // by EPN
  float fiducialVolume;                        // 0.9 - jetR
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  std::vector<fastjet::PseudoJet> leadingJetConstituents;
  JetFinder jetReclusterer;

  HistogramRegistry spectra{
    "spectra",
    {{"vertexZ", "z vertex", {HistType::kTH1F, {{400, -20., +20.}}}}, //
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
     {"fLeadJetEtaVsLeadingTrackEtaPathologicalAll",
      "leading jet eta versus track eta for cases pT jet is less than track pT in TPC volume",
      {HistType::kTH2F, {{40, -1., 1.}, {40, -1., 1.}}}}, //
     {"fLeadJetPhiVsLeadingTrackPhiPathologicalAll",
      "leading jet phi versus track phi for cases pT jet is less than track pT in TPC volume",
      {HistType::kTH2F,
       {{60, 0, TMath::TwoPi()}, {60, 0, TMath::TwoPi()}}}}, //
     {"fLeadJetChPtVsLeadingTrackPathologicalSameDirectionAll",
      "Leading Jet Pt Versus Leading Track Pt  Pathological cases pointing to same direction",
      {HistType::kTH2F, {{200, 0., +200.}, {200, 0., +200.}}}}, //
     {"numberOfJetConstituentsPathologicalSameDirectionAll",
      "number of jet constituents pathological cases pointing to same direction",
      {HistType::kTH1F, {{100, 0., +100.}}}}, //
     {"fLeadJetChPtVsLeadingTrackPathologicalTrackOnBoarderAll",
      "Leading Jet Pt Versus Leading Track Pt  Pathological cases with track on TPC boarder",
      {HistType::kTH2F, {{200, 0., +200.}, {200, 0., +200.}}}}, //
     {"numberOfJetConstituentsPathologicalTrackOnBoarderAll",
      "number of jet constituents pathological cases with track on TPC boarder",
      {HistType::kTH1F, {{100, 0., +100.}}}}, //
     {"fLeadJetChPtVsLeadingTrackPathologicalBulkAll",
      "Leading Jet Pt Versus Leading Track Pt  Pathological cases elsewhere",
      {HistType::kTH2F, {{200, 0., +200.}, {200, 0., +200.}}}}, //
     {"numberOfJetConstituentsPathologicalBulkAll",
      "number of jet constituents pathological cases elsewhere",
      {HistType::kTH1F, {{100, 0., +100.}}}}, //
     {"fLeadJetEtaVsLeadingTrackEtaPathologicalGap",
      "leading jet eta versus track eta for cases pT jet is less than track pT in TPC volume",
      {HistType::kTH2F, {{40, -1., 1.}, {40, -1., 1.}}}}, //
     {"fLeadJetPhiVsLeadingTrackPhiPathologicalGap",
      "leading jet phi versus track phi for cases pT jet is less than track pT in TPC volume",
      {HistType::kTH2F, {{60, 0, TMath::TwoPi()}, {60, 0, TMath::TwoPi()}}}}, //
     {"jetAreaFullVol",
      "area of all jets in full TPC volume",
      {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}}}, //
     {"jetAreaFullVolPathologicalGap",
      "area of all jets in full TPC volume in pathological events",
      {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}}}, //
     {"jetPtminusJetPerp",
      "jetPt minus JetPerp",
      {HistType::kTH1F, {{200, -10., +10.}}}}}};

  // TrackSelection globalTracks;
  void init(o2::framework::InitContext&)
  {

    // globalTracks = getGlobalTrackSelection();
    // globalTracks.SetEtaRange(-cfgTPCVolume, cfgTPCVolume);

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = cfgJetR;
    jetReclusterer.jetEtaMin = -cfgTPCVolume;
    jetReclusterer.jetEtaMax = cfgTPCVolume;

    fiducialVolume = cfgTPCVolume - cfgJetR;
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,
                                    aod::TracksDCA, aod::TrackSelection>;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVertexCut);

  Filter trackEtaFilter = nabs(aod::track::eta) < cfgTPCVolume;
  Filter trackPtFilter = aod::track::pt > cfgJetPtMin;

  void
    process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                    aod::JetFilters>>::iterator const& collision,
            soa::Filtered<TrackCandidates> const& tracks)
  {

    if (collision.hasJetChHighPt() >= bTriggerDecision) {
      jetConstituents.clear();
      jetReclustered.clear();
      leadingJetConstituents.clear();

      float leadingJetPt = -1.0;
      float leadingJetEta = -2.0;
      float leadingJetPhi = -1.0;
      float leadingTrackPt = -1.0;
      float leadingTrackEta = -2.0;
      float leadingTrackPhi = -1.0;
      float leadingJetArea = -1.0;
      float deta = -1.0;
      float dphi = -1.0;

      spectra.fill(HIST("vertexZ"),
                   collision.posZ()); // Inclusive Track Cross TPC Rows

      for (auto& trk : tracks) { //loop over filtered tracks in full TPC volume having pT > 100 MeV

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

        fillConstituents(
          trk,
          jetConstituents); // ./PWGJE/Core/JetFinder.h
                            // recombination scheme is assumed
                            // to be Escheme with pion mass

        if (trk.pt() >
            leadingTrackPt) { // Find leading track pT in full TPC volume
          leadingTrackPt = trk.pt();
          leadingTrackEta = trk.eta();
          leadingTrackPhi = trk.phi();
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

      // Find leading jet pT in full TPC volume
      for (auto& jet : jetReclustered) {
        if (fabs(jet.eta()) < cfgTPCVolume) {

          spectra.fill(HIST("jetPtminusJetPerp"), jet.pt() - jet.perp());

          if (jet.perp() > leadingJetPt) {
            leadingJetPt = jet.perp();
            leadingJetEta = jet.eta();
            leadingJetPhi = jet.phi();
            leadingJetArea = jet.area();
            leadingJetConstituents.clear();
            leadingJetConstituents = jet.constituents();
          }
        }
      }

      if (leadingJetPt > -1.) {
        spectra.fill(HIST("ptphiLeadingJet"), leadingJetPt, leadingJetPhi);
        spectra.fill(HIST("ptetaLeadingJet"), leadingJetPt, leadingJetEta);

        spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt,
                     leadingJetPt); // leading jet pT versus leading track pT
      }

      if (leadingJetPt > -1. && leadingTrackPt > -1.) {
        deta = fabs(leadingTrackEta - leadingJetEta);
        dphi = fabs(TVector2::Phi_mpi_pi(leadingTrackPhi - leadingJetPhi));

        if ((leadingTrackPt - leadingJetPt) > 1e-4) { // pathological case
          spectra.fill(HIST("fLeadJetEtaVsLeadingTrackEtaPathologicalAll"),
                       leadingTrackEta, leadingJetEta);
          spectra.fill(HIST("fLeadJetPhiVsLeadingTrackPhiPathologicalAll"),
                       leadingTrackPhi, leadingJetPhi);

          if (deta < 1e-4 && dphi < 1e-4) {
            spectra.fill(HIST("fLeadJetChPtVsLeadingTrackPathologicalSameDirectionAll"), leadingTrackPt,
                         leadingJetPt); // leading jet pT versus leading track pT

            spectra.fill(HIST("numberOfJetConstituentsPathologicalSameDirectionAll"), leadingJetConstituents.size());
          } else if (fabs(fabs(leadingTrackEta) - 0.9) < 1e-4) {
            spectra.fill(HIST("fLeadJetChPtVsLeadingTrackPathologicalTrackOnBoarderAll"), leadingTrackPt,
                         leadingJetPt); // leading jet pT versus leading track pT

            spectra.fill(HIST("numberOfJetConstituentsPathologicalTrackOnBoarderAll"), leadingJetConstituents.size());
          } else {
            spectra.fill(HIST("fLeadJetChPtVsLeadingTrackPathologicalBulkAll"), leadingTrackPt,
                         leadingJetPt); // leading jet pT versus leading track pT

            spectra.fill(HIST("numberOfJetConstituentsPathologicalBulkAll"), leadingJetConstituents.size());
          }
        }

        if ((leadingTrackPt - leadingJetPt) > 1.) { // pathological case with larger momentum difference
          spectra.fill(HIST("fLeadJetEtaVsLeadingTrackEtaPathologicalGap"),
                       leadingTrackEta, leadingJetEta);
          spectra.fill(HIST("fLeadJetPhiVsLeadingTrackPhiPathologicalGap"),
                       leadingTrackPhi, leadingJetPhi);
          spectra.fill(HIST("jetAreaFullVolPathologicalGap"),
                       leadingJetPt, leadingJetArea);
        }
      }

      //--------------------------------------------------------------

      // Inclusive Jet pT spectrum in Fiducial volume
      for (auto& jet : jetReclustered) {
        if (fabs(jet.eta()) < fiducialVolume) {
          spectra.fill(HIST("ptJetChPtInclFidVol"), jet.perp());
          spectra.fill(HIST("ptphiJetChPtInclFidVol"), jet.perp(), jet.phi());
          spectra.fill(HIST("ptetaJetChPtInclFidVol"), jet.perp(), jet.eta());
        }
        spectra.fill(HIST("ptetaJetChPtInclFullVol"), jet.perp(), jet.eta());
        spectra.fill(HIST("jetAreaFullVol"), jet.perp(), jet.area());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<ChJetTriggerQATask>(
    cfgc, TaskName{"jet-charged-trigger-qa"})};
}
