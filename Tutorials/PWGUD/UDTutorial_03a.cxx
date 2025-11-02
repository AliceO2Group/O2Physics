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
// \brief UD tutorial
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  October 2023

#include "PWGUD/Core/UDHelpers.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDTutorial03a {
  SliceCache cache;

  // a pdg object
  TDatabasePDG* pdg = nullptr;

  // get a DGCutparHolder
  Configurable<int> verbosity{"Verbosity", 0, "Determines level of verbosity"};

  // initialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCMu, aod::McTrackLabels>;
  using TC = TCs::iterator;

  void init(InitContext& context)
  {
    // PDG
    pdg = TDatabasePDG::Instance();

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMCTruth")) {
      registry.add("MC/Stat", "Count generated events; ; Entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
      registry.add("MC/nParts", "Number of McParticles per collision; Number of McParticles; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("MC/genEtaPt", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("MC/genRap", "Generated events; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/genMPt", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("MC/accEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("MC/accRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/accMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
    }

    if (context.mOptions.get<bool>("processReco")) {
      registry.add("Reco/Stat", "Count reconstruted events; ; Entries", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("Reco/nTracks", "Number of reconstructed tracks per collision; Number of reconstructed tracks; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("Reco/nPVContributors", "Number of PV contributors per collision; Number of PV contributors; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("Reco/selEtaPt", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("Reco/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("Reco/selMPt", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
    }
  }

  Preslice<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;

  // retrieve particle mass (GeV/c^2) from TDatabasePDG
  float particleMass(TDatabasePDG* pdg, int pid)
  {
    auto mass = 0.;
    TParticlePDG* pdgparticle = pdg->GetParticle(pid);
    if (pdgparticle != nullptr) {
      mass = pdgparticle->Mass();
    }
    return mass;
  }

  // check if a reconstructed track is a muon candidate
  bool isMuonCandidate_rec(TC track)
  {
    if (abs(track.tpcNSigmaMu()) > 3.) {
      return false;
    }
    return true;
  }

  // check if a generated event is of the type J/Psi -> mu+ + mu- using the MC particle stack
  template <typename MCTrack>
  std::vector<int64_t> getDaughterParts_gen(MCTrack const& parts)
  {
    std::vector<int64_t> selectedParts;

    // in this case we expect the data files to contain events of the type J/Psi -> mu+ + mu-
    if (udhelpers::isSTARLightJPsimumu(parts)) {
      selectedParts.push_back(1);
      selectedParts.push_back(2);
    }
    return selectedParts;
  }

  // retrieve the two muon candidates of a given reconstructed collision
  std::vector<int64_t> getDaughterTracks_rec(CC const& collision, TCs const& tracks)
  {
    // return a vector of track indices
    std::vector<int64_t> emptySelection;
    std::vector<int64_t> selectedTracks;

    // the collision must have exactly 2 PV tracks
    if (collision.numContrib() != 2) {
      return emptySelection;
    }

    // the 2 PV tracks must have opposite charge signs
    // and be muon candidates
    int netCharge = 0;
    int ind = -1;
    for (auto track : tracks) {
      ind++;
      if (track.isPVContributor()) {
        if (!isMuonCandidate_rec(track)) {
          return emptySelection;
        }
        netCharge += track.sign();
        selectedTracks.push_back(ind);
      }
    }
    if (netCharge != 0) {
      return emptySelection;
    }

    // the two PV contributors are both muon candidates
    return selectedTracks;
  }

  // compute the 2-track invariant mass using muon hypothesis
  // bool computeIVM_gen(aod::McParticles const& parts, std::vector<int64_t> partIds, TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
  template <typename TMcPart>
  bool computeIVM_gen(TMcPart const& parts, std::vector<int64_t> partIds, TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
  {
    if (partIds.size() != 2) {
      return false;
    }

    // first particle
    auto part1 = parts.iteratorAt(partIds.at(0));
    auto m1 = particleMass(pdg, 13);
    auto ene1 = sqrt(pow(part1.px(), 2.) + pow(part1.py(), 2.) + pow(part1.pz(), 2.) + pow(m1, 2.));
    *lv1 = TLorentzVector(part1.px(), part1.py(), part1.pz(), ene1);

    // second particle
    auto part2 = parts.iteratorAt(partIds.at(1));
    auto m2 = particleMass(pdg, 13);
    auto ene2 = sqrt(pow(part2.px(), 2.) + pow(part2.py(), 2.) + pow(part2.pz(), 2.) + pow(m2, 2.));
    *lv2 = TLorentzVector(part2.px(), part2.py(), part2.pz(), ene2);

    // system
    *lv = *lv1 + *lv2;

    return true;
  }

  bool computeIVM_rec(TCs const& tracks, std::vector<int64_t> trackIds, TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
  {
    if (trackIds.size() != 2) {
      return false;
    }

    // first track
    auto tr1 = tracks.iteratorAt(trackIds.at(0));
    auto m1 = particleMass(pdg, 13);
    auto ene1 = sqrt(pow(tr1.px(), 2.) + pow(tr1.py(), 2.) + pow(tr1.pz(), 2.) + pow(m1, 2.));
    *lv1 = TLorentzVector(tr1.px(), tr1.py(), tr1.pz(), ene1);

    // second track
    auto tr2 = tracks.iteratorAt(trackIds.at(1));
    auto m2 = particleMass(pdg, 13);
    auto ene2 = sqrt(pow(tr2.px(), 2.) + pow(tr2.py(), 2.) + pow(tr2.pz(), 2.) + pow(m2, 2.));
    *lv2 = TLorentzVector(tr2.px(), tr2.py(), tr2.pz(), ene2);

    // system
    *lv = *lv1 + *lv2;

    return true;
  }

  // check a pair of particles to be in the acceptance of the detector
  bool isInAcceptance(TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
  {
    // the eta of the two muons to be in [-0.9, 0.9]
    // the pT of the two muons to be in [0.1, infty]
    // the y of the (mu+ + mu-)-system to be in [-0.9, 0.9]

    // check first track
    if (abs(lv1->Eta()) > 0.9 || lv1->Pt() < 0.1) {
      return false;
    }

    // check second track
    if (abs(lv2->Eta()) > 0.9 || lv2->Pt() < 0.1) {
      return false;
    }

    // check system
    if (abs(lv->Rapidity()) > 0.9) {
      return false;
    } else {
      return true;
    }
  }

  // check a reconstructed pair of tracks to be a candidate for an event of the type J/Psi -> mu+ + mu-
  bool isSelected_rec(TCs const& tracks, std::vector<int64_t> const& trackIds)
  {
    // tracks is expected to contain two tracks
    if (trackIds.size() != 2) {
      return false;
    }

    // first track
    auto tr1 = tracks.iteratorAt(trackIds.at(0));
    if (!tr1.isPVContributor()) {
      return false;
    }

    // second track
    auto tr2 = tracks.iteratorAt(trackIds.at(1));
    if (!tr2.isPVContributor()) {
      return false;
    }

    // check tracks to be muon candidates
    if (!isMuonCandidate_rec(tr1)) {
      return false;
    }
    if (!isMuonCandidate_rec(tr2)) {
      return false;
    }

    // apply selection criteria
    TLorentzVector* lv1 = new TLorentzVector();
    TLorentzVector* lv2 = new TLorentzVector();
    TLorentzVector* lv = new TLorentzVector();
    if (!computeIVM_rec(tracks, trackIds, lv1, lv2, lv)) {
      return false;
    }

    // select generated events which are in the acceptance
    // (!isInAcceptance(lv1, lv2, lv)) {
    //  return false;
    //}

    // more selection criteria can be applied here ...

    // has passed all selection crtieria
    return true;
  }

  // ...............................................................................................................
  void processMCTruth(aod::McCollisions const& mccollisions, CCs const& /*collisions*/, aod::McParticles const& McParts, TCs const& /*tracks*/)
  {
    // number of McCollisions in DF
    if (verbosity > 0) {
      LOGF(info, "Number of MC collisions %d", mccollisions.size());
    }

    // some variables
    TLorentzVector* lv1_gen = new TLorentzVector();
    TLorentzVector* lv2_gen = new TLorentzVector();
    TLorentzVector* lv_gen = new TLorentzVector();

    // loop over all genererated collisions
    for (auto mccollision : mccollisions) {
      registry.get<TH1>(HIST("MC/Stat"))->Fill(0., 1.);

      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/nParts"))->Fill(partSlice.size(), 1.);
      if (verbosity > 0) {
        LOGF(info, "Number of McParts %d", partSlice.size());
      }

      // compute M versus Pt distributions for 2 cases
      // 1. MC/genMPt - using all generated events, McTruth values
      // 2. MC/accMPt - generated events which are within detector acceptance, McTruth values

      // select generated events of interest
      // retrieve the index of the daughter particles in partSlice
      // we expect there to be exactly 2
      auto daughterPartIds = getDaughterParts_gen(partSlice);
      if (daughterPartIds.size() != 2) {
        continue;
      }

      // use the two selected McParticles to compute the invariant mass
      // lv*_gen are TLorentzVectors of the two tracks and the sum
      if (!computeIVM_gen(partSlice, daughterPartIds, lv1_gen, lv2_gen, lv_gen)) {
        continue;
      }
      registry.get<TH2>(HIST("MC/genEtaPt"))->Fill(lv1_gen->Eta(), lv1_gen->Pt(), 1.);
      registry.get<TH2>(HIST("MC/genEtaPt"))->Fill(lv2_gen->Eta(), lv2_gen->Pt(), 1.);
      registry.get<TH1>(HIST("MC/genRap"))->Fill(lv_gen->Rapidity(), 1.);
      registry.get<TH2>(HIST("MC/genMPt"))->Fill(lv_gen->M(), lv_gen->Pt(), 1.);
      if (verbosity > 0) {
        LOGF(info, "IVMgen %f pT %f", lv_gen->M(), lv_gen->Pt());
      }

      // select generated events which are in the acceptance
      if (!isInAcceptance(lv1_gen, lv2_gen, lv_gen)) {
        continue;
      }
      registry.get<TH2>(HIST("MC/accEtaPt"))->Fill(lv1_gen->Eta(), lv1_gen->Pt(), 1.);
      registry.get<TH2>(HIST("MC/accEtaPt"))->Fill(lv2_gen->Eta(), lv2_gen->Pt(), 1.);
      registry.get<TH1>(HIST("MC/accRap"))->Fill(lv_gen->Rapidity(), 1.);
      registry.get<TH2>(HIST("MC/accMPt"))->Fill(lv_gen->M(), lv_gen->Pt(), 1.);
    }
  }
  PROCESS_SWITCH(UDTutorial03a, processMCTruth, "Process MC truth", true);

  // ...............................................................................................................
  void processReco(CC const& collision, TCs const& tracks, aod::McCollisions const& /*mccollisions*/, aod::McParticles const& /*McParts*/)
  {
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(0., 1.);
    registry.get<TH1>(HIST("Reco/nTracks"))->Fill(tracks.size(), 1.);
    int nContributors = 0;
    for (auto track : tracks) {
      if (track.isPVContributor()) {
        nContributors++;
      }
    }
    registry.get<TH1>(HIST("Reco/nPVContributors"))->Fill(nContributors, 1.);

    // some variables
    TLorentzVector* lv1_rec = new TLorentzVector();
    TLorentzVector* lv2_rec = new TLorentzVector();
    TLorentzVector* lv_rec = new TLorentzVector();

    // select 2 muon candidates using the reconstructed information
    // MC truth is not considered in this step
    auto daughterTrackIds = getDaughterTracks_rec(collision, tracks);
    if (daughterTrackIds.size() != 2) {
      if (verbosity > 0) {
        LOGF(info, "Found %d daughter tracks.", daughterTrackIds.size());
      }
      return;
    }
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(1., 1.);

    // apply all selection cuts on the selected reconstructed tracks
    if (!isSelected_rec(tracks, daughterTrackIds)) {
      return;
    }
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(2., 1.);

    // compute the invariant mass using the reconstructed information
    if (!computeIVM_rec(tracks, daughterTrackIds, lv1_rec, lv2_rec, lv_rec)) {
      return;
    }
    registry.get<TH2>(HIST("Reco/selEtaPt"))->Fill(lv1_rec->Eta(), lv1_rec->Pt(), 1.);
    registry.get<TH2>(HIST("Reco/selEtaPt"))->Fill(lv2_rec->Eta(), lv2_rec->Pt(), 1.);
    registry.get<TH1>(HIST("Reco/selRap"))->Fill(lv_rec->Rapidity(), 1.);
    registry.get<TH2>(HIST("Reco/selMPt"))->Fill(lv_rec->M(), lv_rec->Pt(), 1.);
  }
  PROCESS_SWITCH(UDTutorial03a, processReco, "Process reconstructed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial03a>(cfgc, TaskName{"udtutorial03a"})};
}
