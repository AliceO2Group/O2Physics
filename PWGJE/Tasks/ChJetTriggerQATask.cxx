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
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/TableProducer/jetfinder.h"

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
    "cfgJetPtMin", 0.1,
    "minimum jet pT constituent cut"}; // minimum jet constituent pT

  Configurable<float> cfgTrackPhiMinCut{"cfgTrackPhiMinCut", -999, "track min phi cut"};
  Configurable<float> cfgTrackPhiMaxCut{"cfgTrackPhiMaxCut", 999, "track max phi cut"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<int> bTriggerDecision{
    "bTriggerDecision", 0,
    "Charged Jet Trigger Decision Selection"}; // 0=MB Event, 1=Event selected
                                               // by EPN

  Configurable<float> cfgPtThr{
    "cfgPtThr", 10.,
    "jet pT threshold for some QA plots"}; // jet pT threshold for some QA plots

  float fiducialVolume; // 0.9 - jetR

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
      {"phietaTrackInclGoodAll",
       "phi vs eta all inclusive good tracks",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"phietaTrackInclGoodHighPt",
       "phi vs eta inclusive good tracks with pT > 10 GeV",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptJetChInclFidVol",
       "inclusive charged jet pT in fiducial volume",
       {HistType::kTH1F, {{200, 0., +200.}}}}, //
      {"ptJetChInclFullVol",
       "inclusive charged jet pT in full volume",
       {HistType::kTH1F, {{200, 0., +200.}}}}, //
      {"ptphiJetChInclFidVol",
       "inclusive charged jet pT vs phi in fiducial volume",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptphiJetChInclFullVol",
       "inclusive charged jet pT vs phi in full TPC volume",
       {HistType::kTH2F, {{100, 0., +100.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaJetChInclFidVol",
       "inclusive charged jet pT vs eta in fiducial volume",
       {HistType::kTH2F, {{100, 0., +100.}, {80, -1., 1.}}}}, //
      {"phietaJetChInclFidVol",
       "inclusive charged jet phi vs eta in fiducial volume",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"phietaJetChInclFullVol",
       "inclusive charged jet phi vs eta in full TPC volume",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"phietaJetChInclFidVolHighPt",
       "inclusive charged jet phi vs eta in fiducial volume",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"phietaJetChInclFullVolHighPt",
       "inclusive charged jet phi vs eta in full TPC volume",
       {HistType::kTH2F, {{80, -1., 1.}, {60, 0, TMath::TwoPi()}}}}, //
      {"ptetaJetChInclFullVol",
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
      {"jetAreaFullVol",
       "area of all jets in full TPC volume",
       {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}}}, //
      {"jetAreaFidVol",
       "area of all jets in fiducial volume",
       {HistType::kTH2F, {{100, 0., +100.}, {50, 0., 2.}}}}, //
      {"jetConsituentsThatWereNotFoundAmongTracksPtEtaPhi",
       "PtEtaPhi of jetConsituentsThatWereNotFoundAmongTracks in full volume",
       {HistType::kTH3F, {{100, 0., +100.}, {40, -1., 1.}, {60, 0., TMath::TwoPi()}}}}, //
      {"jetConsituentsThatWereFoundAmongTracksPtEtaPhi",
       "PtEtaPhi of jetConsituentsThatWereFoundAmongTracks in full volume",
       {HistType::kTH3F, {{100, 0., +100.}, {40, -1., 1.}, {60, 0., TMath::TwoPi()}}}}, //
      {"tracksThatWereNotJetConstituentsPtEtaPhi",
       "PtEtaPhi of tracksThatWereNotjetConsituents in full volume",
       {HistType::kTH3F, {{100, 0., +100.}, {40, -1., 1.}, {60, 0., TMath::TwoPi()}}}}, //
      {"tracksThatWereJetConstituentsPtEtaPhi",
       "PtEtaPhi of tracksThatWerejetConsituents in full volume",
       {HistType::kTH3F, {{100, 0., +100.}, {40, -1., 1.}, {60, 0., TMath::TwoPi()}}}} //
    }};

  int eventSelection = -1;
  int trackSelection = -1;
  void init(o2::framework::InitContext&)
  {
    fiducialVolume = static_cast<float>(cfgTPCVolume) - static_cast<float>(cfgJetR);
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(evSel));
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  // declare filters on collisions
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < static_cast<float>(cfgVertexCut));

  // declare filters on tracks
  Filter trackFilter = (nabs(aod::jtrack::eta) < static_cast<float>(cfgTPCVolume)) && (aod::jtrack::phi > static_cast<float>(cfgTrackPhiMinCut)) && (aod::jtrack::phi < static_cast<float>(cfgTrackPhiMaxCut)) && (aod::jtrack::pt > static_cast<float>(cfgJetPtMin));

  // declare filters on jets
  Filter jetRadiusSelection = o2::aod::jet::r == nround(cfgJetR.node() * 100.0f);

  using filteredJets = o2::soa::Filtered<o2::aod::ChargedJets>;
  using TrackCandidates = aod::JTracks;

  void
    process(soa::Filtered<soa::Join<aod::JCollisions,
                                    aod::JChTrigSels>>::iterator const& collision,
            soa::Filtered<TrackCandidates> const& tracks, o2::soa::Filtered<soa::Join<o2::aod::ChargedJets, aod::ChargedJetConstituents>> const& jets)
  // soa::Filtered<TrackCandidates> const& tracks, filteredJets const& jets)
  {

    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }

    if (JetDerivedDataUtilities::selectChargedTrigger(collision, JetDerivedDataUtilities::JTrigSelCh::charged) >= static_cast<int>(bTriggerDecision)) {

      float leadingJetPt = -1.0;
      float leadingJetEta = -2.0;
      float leadingJetPhi = -1.0;
      float leadingTrackPt = -1.0;
      float leadingTrackEta = -2.0;
      float leadingTrackPhi = -1.0;

      spectra.fill(HIST("vertexZ"),
                   collision.posZ()); // Inclusive Track Cross TPC Rows

      std::vector<TLorentzVector> acceptedTracks;
      std::vector<bool> isJetConstituent;
      acceptedTracks.resize(0);
      isJetConstituent.resize(0);

      TLorentzVector v;

      for (auto& trk : tracks) { //loop over filtered tracks in full TPC volume having pT > 100 MeV

        if (!JetDerivedDataUtilities::selectTrack(trk, trackSelection)) {
          continue;
        }
        v.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), 0.139);
        acceptedTracks.push_back(v);
        isJetConstituent.push_back(false); // initialization

        spectra.fill(
          HIST("ptphiTrackInclGood"), trk.pt(),
          trk.phi()); // Inclusive Track pT vs phi spectrum in TPC volume
        spectra.fill(
          HIST("ptetaTrackInclGood"), trk.pt(),
          trk.eta()); // Inclusive Track pT vs eta spectrum in TPC volume

        spectra.fill(
          HIST("phietaTrackInclGoodAll"), trk.eta(),
          trk.phi()); // Inclusive Track pT vs eta spectrum in TPC volume

        if (trk.pt() > static_cast<float>(cfgPtThr)) {
          spectra.fill(
            HIST("phietaTrackInclGoodHighPt"), trk.eta(),
            trk.phi()); // Inclusive Track pT vs eta spectrum in TPC volume
        }

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

      // Find leading jet pT in full TPC volume
      for (auto& jet : jets) {
        if (fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {

          if (jet.pt() > leadingJetPt) {
            leadingJetPt = jet.pt();
            leadingJetEta = jet.eta();
            leadingJetPhi = jet.phi();
          }

          // access jet constituents as tracks
          for (auto& jct : jet.tracks_as<TrackCandidates>()) {
            bool wasFoundAmongTracks = false;
            for (UInt_t itr = 0; itr < acceptedTracks.size(); itr++) {
              if (TMath::Abs(acceptedTracks[itr].Pt() - jct.pt()) < 1e-4 && TMath::Abs(acceptedTracks[itr].Eta() - jct.eta()) < 1e-4 && TMath::Abs(TVector2::Phi_mpi_pi(acceptedTracks[itr].Phi() - jct.phi()) < 1e-4)) {
                isJetConstituent[itr] = true; // initialization
                wasFoundAmongTracks = true;
                break;
              }
            }
            if (!wasFoundAmongTracks) {
              spectra.fill(HIST("jetConsituentsThatWereNotFoundAmongTracksPtEtaPhi"), jct.pt(), jct.eta(), TVector2::Phi_0_2pi(jct.phi()));
            } else {
              spectra.fill(HIST("jetConsituentsThatWereFoundAmongTracksPtEtaPhi"), jct.pt(), jct.eta(), TVector2::Phi_0_2pi(jct.phi()));
            }
          }
        }
      }

      for (UInt_t itr = 0; itr < acceptedTracks.size(); itr++) {
        if (!isJetConstituent[itr]) {
          spectra.fill(HIST("tracksThatWereNotJetConstituentsPtEtaPhi"), acceptedTracks[itr].Pt(), acceptedTracks[itr].Eta(), TVector2::Phi_0_2pi(acceptedTracks[itr].Phi()));
        } else {
          spectra.fill(HIST("tracksThatWereJetConstituentsPtEtaPhi"), acceptedTracks[itr].Pt(), acceptedTracks[itr].Eta(), TVector2::Phi_0_2pi(acceptedTracks[itr].Phi()));
        }
      }

      if (leadingJetPt > -1.) {
        spectra.fill(HIST("ptphiLeadingJet"), leadingJetPt, leadingJetPhi);
        spectra.fill(HIST("ptetaLeadingJet"), leadingJetPt, leadingJetEta);
      }

      if (leadingJetPt > -1. && leadingTrackPt > -1.) {
        spectra.fill(HIST("fLeadJetChPtVsLeadingTrack"), leadingTrackPt,
                     leadingJetPt); // leading jet pT versus leading track pT

        if ((leadingTrackPt - leadingJetPt) > 1e-4) { // pathological case
          spectra.fill(HIST("fLeadJetEtaVsLeadingTrackEtaPathologicalAll"),
                       leadingTrackEta, leadingJetEta);
          spectra.fill(HIST("fLeadJetPhiVsLeadingTrackPhiPathologicalAll"),
                       leadingTrackPhi, leadingJetPhi);
        }
      }

      // Inclusive Jet pT spectrum in Fiducial volume
      for (auto& jet : jets) {
        if (fabs(jet.eta()) < fiducialVolume) {
          spectra.fill(HIST("ptJetChInclFidVol"), jet.pt());
          spectra.fill(HIST("ptphiJetChInclFidVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFidVol"), jet.pt(), jet.eta());
          spectra.fill(HIST("phietaJetChInclFidVol"), jet.eta(), jet.phi());
          if (jet.pt() > static_cast<float>(cfgPtThr)) {
            spectra.fill(HIST("phietaJetChInclFidVolHighPt"), jet.eta(), jet.phi());
          }
          spectra.fill(HIST("jetAreaFidVol"), jet.pt(), jet.area());
        }

        if (fabs(jet.eta()) < static_cast<float>(cfgTPCVolume)) {
          spectra.fill(HIST("ptJetChInclFullVol"), jet.pt());
          spectra.fill(HIST("ptphiJetChInclFullVol"), jet.pt(), jet.phi());
          spectra.fill(HIST("ptetaJetChInclFullVol"), jet.pt(), jet.eta());
          spectra.fill(HIST("phietaJetChInclFullVol"), jet.eta(), jet.phi());
          if (jet.pt() > static_cast<float>(cfgPtThr)) {
            spectra.fill(HIST("phietaJetChInclFullVolHighPt"), jet.eta(), jet.phi());
          }
          spectra.fill(HIST("jetAreaFullVol"), jet.pt(), jet.area());
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
