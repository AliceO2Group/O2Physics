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
/// \file jetDsSpecSubs.cxx
/// \brief Ds-tagged jet analysis with substructure histogram outputs
/// \author Monalisa Melo <monalisa.melo@cern.ch>, Universidade de São Paulo

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include "TVector3.h"


#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


consteval float getValFromBin(int bin)
{
  return static_cast<float>(bin) - 0.5f;
}

enum BinExpColCntr {AllCollisions = 1,
                    Sel8ZCut = 2};

enum BinExpJetCntr {ChargedJets = 1 };
enum BinMCColCntr {All = 1,
                   ZCut = 2,
                   Matched = 3,
                   MatchedSel8ZCut = 4
                  };

enum BinMCJetCntr {DetectorLevelJetInMCCollision = 1,
                   ParticleLevelJetInMCCollision = 2,
                   DetectorLevelJetWithMatchedCandidate = 3,
                   ParticleLevelJetWithMatchedCandidate = 4
};



struct JetDsSpecSubs {

  //==================
  // Type definitions
  //==================

  using DsCandidatesData = aod::CandidatesDsData;
  using DsCandidatesMCD = aod::CandidatesDsMCD;
  using DsCandidatesMCP = aod::CandidatesDsMCP;

  using DsDataJets = soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents>;
  using DsMCDJets = soa::Join<aod::DsChargedMCDetectorLevelJets, aod::DsChargedMCDetectorLevelJetConstituents, aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>;
  using DsMCPJets = soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>;

  // Slices for access to proper HF MCD jet collision that is associated to MCCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<DsMCDJets> dsMCDJetsPerEXPCollisionPreslice = aod::jet::collisionId;
  Preslice<DsMCPJets> dsMCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // internals
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  // Filters
  //Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  //Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  //=============
  // Histograms
  //=============

  HistogramRegistry registry{
    "registry",
    { {"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}},
      {"h_collision_counter_data", ";event counter;entries", {HistType::kTH1F, {{10, 0., 10.}}}},
      {"h_collision_counter_mcd", ";event counter;entries", {HistType::kTH1F, {{10, 0., 10.}}}},
      {"h_collision_counter_mcp", ";event counter;entries", {HistType::kTH1F, {{10, 0., 10.}}}},

      {"h_track_pt", ";#it{p}_{T,track};entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

      // Data histograms
      {"h_dsjet_counter_data", ";type;counts", {HistType::kTH1F, {{3, 0., 3.}}}},

      {"h_jet_pt_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_data", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_data", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_mass_data", ";m_{D_{S}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{1000, 0., 6.}}}},
      {"h_ds_pt_data", ";#it{p}_{T,D_{S}} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_eta_data", ";#eta_{D_{S}};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_phi_data", ";#phi_{D_{S}};entries", {HistType::kTH1F, {{250, -1., 7.}}}},

      {"h_ds_jet_pt_data", ";#it{p}_{T,D_{S} jet} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_jet_eta_data", ";#eta_{D_{S} jet};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_jet_phi_data", ";#phi_{D_{S} jet};entries", {HistType::kTH1F, {{250, -1., 7.}}}},

      {"h_ds_jet_projection_data", ";z^{D_{S},jet}_{||};entries", {HistType::kTH1F, {{1000, 0., 2.}}}},
      {"h_ds_jet_distance_data", ";#DeltaR_{D_{S},jet};entries", {HistType::kTH1F, {{1000, 0., 1.}}}},
      {"h_ds_jet_mass_data", ";m_{jet}^{ch} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
      {"h_ds_jet_lambda11_data", ";#lambda_{1}^{1};entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"h_ds_jet_lambda12_data", ";#lambda_{2}^{1};entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"hSparse_ds_data", ";m_{D_{S}};#it{p}_{T,D_{S}};#it{p}_{T,jet};z^{D_{S},jet}_{||};#DeltaR_{D_{S},jet}", {HistType::kTHnSparseF, {{60, 1.7, 2.1}, {60, 0., 100.}, {60, 0., 100.}, {60, 0., 2.}, {60, 0., 1.0}}}},

      // MC detector-level histograms
      {"McEffCol", "N_{collisions};", {HistType::kTH1F, {{4, 0., 4.0}}}},
      {"McEffJet", "N_{jet};", {HistType::kTH1F, {{4, 0., 4.0}}}},
      {"h_jet_pt_mcd", "detector-level jet pT;#it{p}_{T,jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_mcd", "detector-level jet #eta;#eta_{jet}^{det};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_mcd", "detector-level jet #phi;#phi_{jet}^{det};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_mass_mcd", ";m_{D_{S}}^{rec} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{1000, 0., 6.}}}},
      {"h_ds_pt_mcd", ";#it{p}_{T,D_{S}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_eta_mcd", ";#eta_{D_{S}}^{det};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_phi_mcd", ";#phi_{D_{S}}^{det};entries", {HistType::kTH1F, {{250, -1., 7.}}}},

      {"h_ds_jet_pt_mcd", ";#it{p}_{T,D_{S} jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_jet_eta_mcd", ";#eta_{D_{S} jet}^{det};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_jet_phi_mcd", ";#phi_{D_{S} jet}^{det};entries", {HistType::kTH1F, {{250, -1., 7.}}}},
      {"h_ds_jet_projection_mcd", ";z^{D_{S},jet}_{||, det};entries", {HistType::kTH1F, {{1000, 0., 2.}}}},
      {"h_ds_jet_distance_mcd", ";#DeltaR_{D_{S},jet}^{det};entries", {HistType::kTH1F, {{1000, 0., 1.}}}},
      {"h_ds_jet_mass_mcd", ";m_{jet}^{ch, det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
      {"h_ds_jet_lambda11_mcd", ";#lambda_{1}^{1, det};entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"h_ds_jet_lambda12_mcd", ";#lambda_{2}^{1, det};entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"hSparse_ds_mcd", ";m_{D_{S}}^{rec};#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||, det};#DeltaR_{D_{S},jet}^{det}", {HistType::kTHnSparseF, {{60, 1.7, 2.1}, {60, 0., 100.}, {60, 0., 100.}, {60, 0., 2.}, {60, 0., 1.0}}}},

      // MC particle-level histograms
      {"h_dsjet_counter_mcp", ";type;counts", {HistType::kTH1F, {{3, 0., 3.}}}},
      {"h_jet_pt_mcp", "particle-level jet pT;#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_mcp", "particle-level jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_mcp", "particle-level jet #phi;#phi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_pt_mcp", ";#it{p}_{T,D_{S}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_eta_mcp", ";#eta_{D_{S}}^{part};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_phi_mcp", ";#phi_{D_{S}}^{part};entries", {HistType::kTH1F, {{250, -1., 7.}}}},

      {"h_ds_jet_pt_mcp", ";#it{p}_{T,D_{S} jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_jet_eta_mcp", ";#eta_{D_{S} jet}^{part};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_jet_phi_mcp", ";#phi_{D_{S} jet}^{part};entries", {HistType::kTH1F, {{250, -1., 7.}}}},
      {"h_ds_jet_projection_mcp", ";z^{D_{S},jet}_{||, part};entries", {HistType::kTH1F, {{1000, 0., 2.}}}},
      {"h_ds_jet_distance_mcp", ";#DeltaR_{D_{S},jet}^{part};entries", {HistType::kTH1F, {{1000, 0., 1.}}}},
      {"hSparse_ds_mcp", ";#it{p}_{T,D_{S}}^{part};#it{p}_{T,jet}^{part};z^{D_{S},jet}_{||, part};#DeltaR_{D_{S},jet}^{part}", {HistType::kTHnSparseF, {{60, 0., 100.}, {60, 0., 100.}, {60, 0., 2.}, {60, 0., 1.0}}}},
    }};
  //========
  // INIT
  //========

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto hData = registry.get<TH1>(HIST("h_dsjet_counter_data"));
    hData->GetXaxis()->SetBinLabel(1, "Ds-jet entries");
    hData->GetXaxis()->SetBinLabel(2, "Ds candidates");
    hData->GetXaxis()->SetBinLabel(3, "Ds jets with >=1 cand.");

    auto hMcp = registry.get<TH1>(HIST("h_dsjet_counter_mcp"));
    hMcp->GetXaxis()->SetBinLabel(1, "Ds-jet entries");
    hMcp->GetXaxis()->SetBinLabel(2, "Ds particles");
    hMcp->GetXaxis()->SetBinLabel(3, "Ds jets with >=1 particle");

     // Labels
    auto mcCollisionCounter = registry.get<TH1>(HIST("McEffCol"));
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::All, "mccollisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::ZCut, "z_cut");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::Matched, "collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::MatchedSel8ZCut, "sel8");

    auto jetCounter = registry.get<TH1>(HIST("McEffJet"));
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetInMCCollision, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetInMCCollision, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetWithMatchedCandidate, "particle matched jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetWithMatchedCandidate, "detector matched jets");

  }
  //===============
  // Lambda compute
  //===============

  template <typename JET, typename TRACKS>
  float computeLambda(JET const& jet, TRACKS const& tracks, float alpha, float kappa)
  {
    if (jet.pt() <= 0.f) {
      return -1.f;
    }

    float sum = 0.f;
    for (auto const& trk : tracks) {
      const float dr = jetutilities::deltaR(jet, trk);
      sum += std::pow(trk.pt(), kappa) * std::pow(dr, alpha);
    }
    const float jetR = jet.r() / 100.f;

    const float denom = std::pow(jet.pt(), kappa) * std::pow(jetR, alpha);
    if (denom <= 0.f) {
      return -1.f;
    }
    return sum / denom;
  }

  //=================
  // Jet Mass compute
  //=================

  template <typename TRACKS>
  float computeJetMass(TRACKS const& tracks)
  {
    double sumPx = 0.0, sumPy = 0.0, sumPz = 0.0, sumE = 0.0;

    for (auto const& trk : tracks) {
      const double pt = trk.pt();
      const double phi = trk.phi();
      const double eta = trk.eta();

      const double px = pt * std::cos(phi);
      const double py = pt * std::sin(phi);
      const double pz = pt * std::sinh(eta);
      const double p = std::sqrt(px * px + py * py + pz * pz);

      sumPx += px;
      sumPy += py;
      sumPz += pz;
      sumE += p; // massless
    }

    const double m2 = sumE * sumE - (sumPx * sumPx + sumPy * sumPy + sumPz * sumPz);
    return (m2 > 0.0) ? static_cast<float>(std::sqrt(m2)) : 0.f;
  }

  //==============
  // Collision QA
  //==============

  void processCollisions(aod::JetCollision const& collision, aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_collisions"), 1.5);

    for (auto const& track : tracks) {

      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        return;
      }

      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }
  PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "collision QA", false);


  //==============
  // DATA process
  //==============

  void processDataChargedSubstructure(aod::JetCollision const& collision, soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents> const& jets,
                                      aod::CandidatesDsData const&, aod::JetTracks const&)
  {
    registry.fill(HIST("h_collision_counter_data"), 2.0);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
        !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }

    registry.fill(HIST("h_collision_counter_data"), 3.0);

    for (const auto& jet : jets) {
      registry.fill(HIST("h_dsjet_counter_data"), 0.5); //DsChargedJets entries

      registry.fill(HIST("h_jet_pt_data"), jet.pt());
      registry.fill(HIST("h_jet_eta_data"), jet.eta());
      registry.fill(HIST("h_jet_phi_data"), jet.phi());

      auto jetTracks = jet.tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);
      const float lambda12 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);


      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      int nDsInJet = 0;

      // Loop over Ds candidates (particle level)
      for (const auto& dsCandidate : jet.candidates_as<aod::CandidatesDsData>()) {
        ++nDsInJet;
        registry.fill(HIST("h_dsjet_counter_data"), 1.5); // Ds candidates associated with the jet

        TVector3 dsVector(dsCandidate.px(), dsCandidate.py(), dsCandidate.pz());

        // Axis distance Delta_R
        const float deltaR = jetutilities::deltaR(jet, dsCandidate);
        // zParallel defined as longitudinal momentum fraction along the jet axis
        const float zParallel = (jetVector * dsVector) / (jetVector * jetVector);

        // --- Ds-level observables ---
        registry.fill(HIST("h_ds_mass_data"), dsCandidate.m());
        registry.fill(HIST("h_ds_pt_data"), dsCandidate.pt());
        registry.fill(HIST("h_ds_eta_data"), dsCandidate.eta());
        registry.fill(HIST("h_ds_phi_data"), dsCandidate.phi());

        registry.fill(HIST("h_ds_jet_projection_data"), zParallel);

        registry.fill(HIST("h_ds_jet_distance_data"), deltaR);

        // Main THnSparse: invariant mass, pT, z, and DeltaR
        registry.fill(HIST("hSparse_ds_data"),
                      dsCandidate.m(),
                      dsCandidate.pt(),
                      jet.pt(),
                      zParallel,
                      deltaR);
      }

      // Jet-level quantities (filled once per jet containing at least one Ds)
      if (nDsInJet > 0) {

        registry.fill(HIST("h_dsjet_counter_data"), 2.5); // Ds jets with at least one associated candidate

        // Jet properties
        registry.fill(HIST("h_ds_jet_pt_data"), jet.pt());
        registry.fill(HIST("h_ds_jet_eta_data"), jet.eta());
        registry.fill(HIST("h_ds_jet_phi_data"), jet.phi());
        // Jet Mass
        registry.fill(HIST("h_ds_jet_mass_data"), mjet);

        // Jet substructure observables
        if (lambda11 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda11_data"), lambda11);
        }
        if (lambda12 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda12_data"), lambda12);
        }

      }
    }
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "Data charged jets", false);


  //==============
  //  MC function
  //==============
  template <typename MCDJetsPerMCCollissionPreslice,
            typename DsMCDJets,
            typename CandidatesMCD,
            typename CandidatesMCP>
  void analyseMonteCarloEfficiency(MCDJetsPerMCCollissionPreslice const& jetmcdpreslice,
                                   aod::JetMcCollisions const& mccollisions,
                                   aod::JetCollisionsMCD const& collisions,
                                   DsMCDJets const& mcdjets,
                                   CandidatesMCD const& /*mcdCandidates*/,
                                   CandidatesMCP const& /*mcpCandidates*/,
                                   aod::JetTracks const& tracks)
  {
    for (const auto& mccollision : mccollisions) {
      // Count all generated MC collisions
      registry.fill(HIST("McEffCol"), getValFromBin(BinMCColCntr::All));

      // Apply MC vertex selection
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      // MC collisions passing z_cut selection
      registry.fill(HIST("McEffCol"), getValFromBin(BinMCColCntr::ZCut));

      // Reconstructed collisions associated to this mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        // Successfully matched reconstructed collision
        registry.fill(HIST("McEffCol"), getValFromBin(BinMCColCntr::Matched));

        // Apply standard event selection and vertex cut
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
            !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        // Matched collision passing analysis selections
        registry.fill(HIST("McEffCol"), getValFromBin(BinMCColCntr::MatchedSel8ZCut));

        // Detector-level Ds-tagged jets associated with the current reconstructed collision
        const auto dsmcdJetsPerCollision = mcdjets.sliceBy(jetmcdpreslice, collision.globalIndex());
        for (const auto& mcdjet : dsmcdJetsPerCollision) {

          // Detector-level jet found in a matched collision
          registry.fill(HIST("McEffJet"), getValFromBin(BinMCJetCntr::DetectorLevelJetInMCCollision));

          // Leading Ds candidate associated to the jet
          auto mcdDscand = mcdjet.template candidates_first_as<CandidatesMCD>();

          // Check whether a matched particle-level jet exists
          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("McEffJet"), getValFromBin(BinMCJetCntr::DetectorLevelJetWithMatchedCandidate));
          }

          // Compute jet-substructure observables
          TVector3 mcd_jetvector(mcdjet.px(), mcdjet.py(), mcdjet.pz());
          TVector3 mcd_candvector(mcdDscand.px(), mcdDscand.py(), mcdDscand.pz());

          float mcd_zParallel = (mcd_jetvector * mcd_candvector) / (mcd_jetvector * mcd_jetvector);
          // Axis distance Delta_R
          float mcd_deltaR = jetutilities::deltaR(mcdjet, mcdDscand);
          float mcd_lambda11 = computeLambda(mcdjet, tracks, 1.f, 1.f);
          float mcd_lambda12 = computeLambda(mcdjet, tracks, 2.f, 1.f);

          // Detector-level Jet Histograms
          registry.fill(HIST("h_jet_pt_mcd"), mcdjet.pt());
          registry.fill(HIST("h_jet_eta_mcd"), mcdjet.eta());
          registry.fill(HIST("h_jet_phi_mcd"), mcdjet.phi());
          registry.fill(HIST("h_ds_jet_projection_mcd"), mcd_zParallel);
          registry.fill(HIST("h_ds_jet_lambda11_mcd"), mcd_lambda11);
          registry.fill(HIST("h_ds_jet_lambda12_mcd"), mcd_lambda12);
          // Detector-level Ds Histgrams
          registry.fill(HIST("h_ds_jet_pt_mcd"), mcdDscand.pt());
          registry.fill(HIST("h_ds_jet_mass_mcd"), mcdDscand.m());
          registry.fill(HIST("h_ds_jet_eta_mcd"), mcdDscand.eta());
          registry.fill(HIST("h_ds_jet_phi_mcd"), mcdDscand.phi());

          // Main THnSparse: invariant mass, pT, z, and DeltaR
          registry.fill(HIST("hSparse_ds_mcd"),
                        mcdDscand.m(),
                        mcdDscand.pt(),
                        mcdjet.pt(),
                        mcd_zParallel,
                        mcd_deltaR);

        }

      }
    }
  }
  //==============
  // MC process
  //==============

  void processMonteCarloEfficiencyDs(aod::JetMcCollisions const& mccollisions,
                                     aod::JetCollisionsMCD const& collisions,
                                     DsMCDJets const& mcdjets,
                                     DsMCPJets const& mcpjets,
                                     aod::CandidatesDsMCD const& mcdCandidates,
                                     aod::CandidatesDsMCP const& mcpCandidates,
                                     aod::JetTracks const& jettracks)
  {

    analyseMonteCarloEfficiency<Preslice<DsMCDJets>,
                                DsMCDJets,
                                DsMCPJets,
                                aod::CandidatesDsMCD,
                                aod::CandidatesDsMCP>(dsMCDJetsPerEXPCollisionPreslice,
                                                      mccollisions,
                                                      collisions,
                                                      mcdjets,
                                                      mcpjets,
                                                      mcdCandidates,
                                                      mcpCandidates,
                                                      jettracks);


  }
  PROCESS_SWITCH(JetDsSpecSubs, processMonteCarloEfficiencyDs, "Non-matched and matched MC Ds and jets", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc)};
}
