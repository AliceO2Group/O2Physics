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

struct UDTutorial03b {
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
      registry.add("MC/recCols", "Number of reconstructed collisions; Number of reconstructed collisions; Entries", {HistType::kTH1F, {{31, -0.5, 30.5}}});
      registry.add("MC/nParts", "Number of McParticles per collision; Number of McParticles; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("MC/nRecTracks", "Number of reconstructed tracks per McParticle; Number of reconstructed tracks per McParticle; Entries", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("MC/genEtaPt", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("MC/genRap", "Generated events; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/genMPt", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("MC/accEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("MC/accRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/accMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("MC/selEtaPt", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("MC/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("MC/selMPt", "Selected events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("MC/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});
    }

    if (context.mOptions.get<bool>("processReco")) {
      registry.add("Reco/Stat", "Count reconstruted events; ; Entries", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("Reco/nTracks", "Number of reconstructed tracks per collision; Number of reconstructed tracks; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("Reco/nPVContributors", "Number of PV contributors per collision; Number of PV contributors; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
      registry.add("Reco/selEtaPt", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("Reco/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("Reco/selMPt", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("Reco/mcEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
      registry.add("Reco/mcRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{300, -1.5, 1.5}}});
      registry.add("Reco/mcMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
      registry.add("Reco/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});
    }
  }

  Preslice<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CCs> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::mctracklabel::mcParticleId;

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

  // check if a reconstructed track represents a muon candidate
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

  // find the McParticles belongin to given tracks
  template <typename MCTrack>
  std::vector<int64_t> getDaughterParts_rec(TCs const& tracks, std::vector<int64_t> trackIds, MCTrack const& /*parts*/)
  {
    std::vector<int64_t> emptySelection;
    std::vector<int64_t> selectedParts;

    for (auto trackId : trackIds) {
      auto tr = tracks.iteratorAt(trackId);
      if (!tr.has_mcParticle()) {
        return emptySelection;
      } else {
        selectedParts.push_back(tr.mcParticle().globalIndex());
      }
    }
    return selectedParts;
  }

  // retrieve the reconstructed tracks which are associated with the given McParticles
  template <typename McPart>
  std::vector<int64_t> getDaughterTracks_gen(McPart const& parts, std::vector<int64_t> partIds, TCs const& tracks)
  {
    // return a vector of track indices
    std::vector<int64_t> emptySelection;
    std::vector<int64_t> selectedTracks;

    for (auto partId : partIds) {
      auto part = parts.iteratorAt(partId);
      auto trs = tracks.sliceBy(trackPerMcParticle, part.globalIndex());
      if (trs.size() == 0) {
        return emptySelection;
      }
      if (trs.size() > 1) {
        LOGF(info, "%d tracks belong to same McParticle!", trs.size());
      }
      for (auto tr : trs) {
        selectedTracks.push_back(tr.globalIndex());
      }
    }
    return selectedTracks;
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
  template <typename TTrack>
  bool computeIVM(TTrack const& tracks, std::vector<int64_t> trackIds, TLorentzVector* lv1, TLorentzVector* lv2, TLorentzVector* lv)
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

    // apply selection criteria
    TLorentzVector* lv1 = new TLorentzVector();
    TLorentzVector* lv2 = new TLorentzVector();
    TLorentzVector* lv = new TLorentzVector();
    if (!computeIVM(tracks, trackIds, lv1, lv2, lv)) {
      return false;
    }
    // select generated events which are in the acceptance
    // (!isInAcceptance(lv1, lv2, lv)) {
    //  return false;
    //}

    // check tracks to be muon candidates
    if (!isMuonCandidate_rec(tr1)) {
      return false;
    }
    if (!isMuonCandidate_rec(tr2)) {
      return false;
    }

    // more selection criteria can be applied here ...

    // has passed all selection crtieria
    return true;
  }

  // ...............................................................................................................
  void processMCTruth(aod::McCollisions const& mccollisions, CCs const& collisions, aod::McParticles const& McParts, TCs const& tracks)
  {
    // number of McCollisions in DF
    if (verbosity > 0) {
      LOGF(info, "Number of MC collisions %d", mccollisions.size());
    }

    // some variables
    TLorentzVector* lv1_gen = new TLorentzVector();
    TLorentzVector* lv2_gen = new TLorentzVector();
    TLorentzVector* lv_gen = new TLorentzVector();
    TLorentzVector* lv1_rec = new TLorentzVector();
    TLorentzVector* lv2_rec = new TLorentzVector();
    TLorentzVector* lv_rec = new TLorentzVector();

    // loop over all generated collisions
    for (auto mccollision : mccollisions) {
      registry.get<TH1>(HIST("MC/Stat"))->Fill(0., 1.);

      // get reconstructed collision which belongs to mccollision
      auto colSlice = collisions.sliceBy(colPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/recCols"))->Fill(colSlice.size(), 1.);

      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/nParts"))->Fill(partSlice.size(), 1.);
      if (verbosity > 0) {
        LOGF(info, "Number of McParts %d", partSlice.size());
      }

      // compute M versus Pt distributions for 3 cases
      // 1. MC/genMPt - using all generated events, McTruth values
      // 2. MC/accMPt - generated events which are within detector acceptance, McTruth values
      // 3. MC/selMPt - events which path selection criteria, reconstructed values

      // select generated events of interest
      // retrieve the index of the daughter particles in partSlice
      // we expect there to be exactly 2
      auto daughterPartIds = getDaughterParts_gen(partSlice);
      if (daughterPartIds.size() != 2) {
        continue;
      }

      // use the two selected McParticles to compute the invariant mass
      // lv*_gen are TLorentzVectors of the two tracks and the sum
      if (!computeIVM(partSlice, daughterPartIds, lv1_gen, lv2_gen, lv_gen)) {
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

      // now obtain the reconstructed tracks
      // which are the tracks associated with the McParticles
      // we expect there to be exactly 2
      auto daughterTrackIds = getDaughterTracks_gen(partSlice, daughterPartIds, tracks);
      if (daughterTrackIds.size() != 2) {
        continue;
      }

      // apply all selection cuts on the selected reconstructed tracks
      if (!isSelected_rec(tracks, daughterTrackIds)) {
        continue;
      }

      // compute the invariant mass using the reconstructed tracks
      if (!computeIVM(tracks, daughterTrackIds, lv1_rec, lv2_rec, lv_rec)) {
        continue;
      }
      registry.get<TH2>(HIST("MC/selEtaPt"))->Fill(lv1_rec->Eta(), lv1_rec->Pt(), 1.);
      registry.get<TH2>(HIST("MC/selEtaPt"))->Fill(lv2_rec->Eta(), lv2_rec->Pt(), 1.);
      registry.get<TH1>(HIST("MC/selRap"))->Fill(lv_rec->Rapidity(), 1.);
      registry.get<TH2>(HIST("MC/selMPt"))->Fill(lv_rec->M(), lv_rec->Pt(), 1.);

      // compute the difference between generated and reconstructed particle momentum
      for (auto McPart : partSlice) {
        // get track which corresponds to McPart
        auto trackSlice = tracks.sliceBy(trackPerMcParticle, McPart.globalIndex());
        registry.get<TH1>(HIST("MC/nRecTracks"))->Fill(trackSlice.size(), 1.);

        // are there reconstructed tracks?
        if (trackSlice.size() > 0) {
          for (auto track : trackSlice) {
            auto pTrack = track.p();
            auto pPart = McPart.p();
            auto pDiff = pTrack - pPart;
            registry.get<TH2>(HIST("MC/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
            if (verbosity > 0) {
              LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: %d dP: %f", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess(), track.isPVContributor(), pDiff);
            }
          }
        } else {
          registry.get<TH2>(HIST("MC/pDiff"))->Fill(-5.9, -1, 1.);
          if (verbosity > 0) {
            LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: No dP: nan", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess());
          }
        }
      }
      if (verbosity > 0) {
        LOGF(info, "");
      }
    }
  }
  PROCESS_SWITCH(UDTutorial03b, processMCTruth, "Process MC truth", true);

  // ...............................................................................................................
  void processReco(CC const& collision, TCs const& tracks, aod::McCollisions const& /*mccollisions*/, aod::McParticles const& McParts)
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
    TLorentzVector* lv1_gen = new TLorentzVector();
    TLorentzVector* lv2_gen = new TLorentzVector();
    TLorentzVector* lv_gen = new TLorentzVector();

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
    if (!computeIVM(tracks, daughterTrackIds, lv1_rec, lv2_rec, lv_rec)) {
      return;
    }
    registry.get<TH2>(HIST("Reco/selEtaPt"))->Fill(lv1_rec->Eta(), lv1_rec->Pt(), 1.);
    registry.get<TH2>(HIST("Reco/selEtaPt"))->Fill(lv2_rec->Eta(), lv2_rec->Pt(), 1.);
    registry.get<TH1>(HIST("Reco/selRap"))->Fill(lv_rec->Rapidity(), 1.);
    registry.get<TH2>(HIST("Reco/selMPt"))->Fill(lv_rec->M(), lv_rec->Pt(), 1.);

    // now access the McTruth information
    // get McCollision belonging to collision
    if (collision.has_mcCollision()) {
      // auto mccollision = collision.mcCollision();
      registry.get<TH1>(HIST("Reco/Stat"))->Fill(3., 1.);
    } else {
      if (verbosity > 0) {
        LOGF(info, "This collision has no associated McCollision");
      }
    }

    // compute the difference between generated and reconstructed momentum
    for (auto track : tracks) {
      // is there an associated McParticle?
      if (track.has_mcParticle()) {
        auto pTrack = track.p();
        auto McPart = track.mcParticle();
        auto pPart = McPart.p();
        auto pDiff = pTrack - pPart;
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
        if (verbosity > 0) {
          LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: %d dP: %f", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess(), track.isPVContributor(), pDiff);
        }
      } else {
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(-5.9, -1, 1.);
        if (verbosity > 0) {
          LOGF(info, "  This track has no associated McParticle");
        }
      }
    }
    if (verbosity > 0) {
      LOGF(info, "");
    }

    // both muon candidates are required to have an associated McParticle
    auto daughterPartIds = getDaughterParts_rec(tracks, daughterTrackIds, McParts);
    if (daughterPartIds.size() != 2) {
      if (verbosity > 0) {
        LOGF(info, "Found %d daughter McParticles.", daughterPartIds.size());
      }
      return;
    } else {
      registry.get<TH1>(HIST("Reco/Stat"))->Fill(4., 1.);

      // compute invariant mass using McTruth
      if (!computeIVM(McParts, daughterPartIds, lv1_gen, lv2_gen, lv_gen)) {
        return;
      }
      registry.get<TH2>(HIST("Reco/mcEtaPt"))->Fill(lv1_gen->Eta(), lv1_gen->Pt(), 1.);
      registry.get<TH2>(HIST("Reco/mcEtaPt"))->Fill(lv2_gen->Eta(), lv2_gen->Pt(), 1.);
      registry.get<TH1>(HIST("Reco/mcRap"))->Fill(lv_gen->Rapidity(), 1.);
      registry.get<TH2>(HIST("Reco/mcMPt"))->Fill(lv_gen->M(), lv_gen->Pt(), 1.);
    }
  }
  PROCESS_SWITCH(UDTutorial03b, processReco, "Process reconstructed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDTutorial03b>(cfgc, TaskName{"udtutorial03b"})};
}
