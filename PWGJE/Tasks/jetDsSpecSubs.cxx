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

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>
#include <TVector3.h>

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

enum BinExpColCntr { AllCollisions = 1,
                     Sel8ZCut = 2 };

enum BinExpJetCntr { ChargedJets = 1 };
enum BinMCColCntr { All = 1,
                    ZCut = 2,
                    Matched = 3,
                    MatchedSel8ZCut = 4
};

enum BinMCJetCntr { DetectorLevelJetInMCCollision = 1,
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
  using DsMCPJets = soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets>;

  using DsDataJetsEWS = soa::Join<
    aod::DsChargedEventWiseSubtractedJets,
    aod::DsChargedEventWiseSubtractedJetConstituents>;

  // Inclusive charged jets
  using ChargedJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
  using ChargedJetsEWS = soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>;

  // Slices for access to proper HF MCD jet collision that is associated to MCCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;

  Preslice<DsMCDJets> dsMCDJetsPerEXPCollisionPreslice = aod::jet::collisionId;
  // Preslice<DsMCDJetsEWS> dsMCDJetsEWSPerEXPCollisionPreslice = aod::jet::collisionId;

  Preslice<DsMCPJets> dsMCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;
  // Preslice<DsMCPJetsEWS> dsMCPJetsEWSPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Event configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // Event-wise constituent subtraction jet tables
  Configurable<float> centralityMin{"centralityMin", -999.f, "Minimum FT0M centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.f, "Maximum FT0M centrality"};

  // internals
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  // Filters
  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  // Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  // Filtered jet tables
  using FilteredDsDataJets = soa::Filtered<DsDataJets>;
  using FilteredDsMCDJets = soa::Filtered<DsMCDJets>;
  using FilteredDsMCPJets = soa::Filtered<DsMCPJets>;

  using FilteredDsDataJetsEWS = soa::Filtered<DsDataJetsEWS>;
  // using FilteredDsMCDJetsEWS = soa::Filtered<DsMCDJetsEWS>;
  // using FilteredDsMCPJetsEWS = soa::Filtered<DsMCPJetsEWS>;

  // Filtered inclusive charged jets
  using FilteredChargedJets = soa::Filtered<ChargedJets>;
  using FilteredChargedJetsEWS = soa::Filtered<ChargedJetsEWS>;

  Filter collisionFilter =
    nabs(aod::jcollision::posZ) < vertexZCut &&
    aod::jcollision::centFT0M >= centralityMin &&
    aod::jcollision::centFT0M < centralityMax;

  //=============
  // Histograms
  //=============

  HistogramRegistry registry{
    "registry",
    {
      {"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}},
      {"h_event_counter_data", ";Selection step;Events", {HistType::kTH1F, {{3, 0.5, 3.5}}}},
      // Track QA
      {"h_track_pt", ";#it{p}_{T,track};entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}}},
      // Inclusive jet QA
      {"h_njets_inclusive", ";Number of charged jets per event;Events", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"h_jet_pt_inclusive", ";#it{p}_{T,jet} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_inclusive", ";#eta_{jet};Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_jet_phi_inclusive", ";#varphi_{jet} (rad);Entries", {HistType::kTH1F, {{80, -1., 7.}}}},
      {"h_jet_mass_inclusive", ";#it{m}_{jet} (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {{120, 0., 60.}}}},
      {"h_jet_nconst_inclusive", ";Jet constituents;Entries", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"h_lambda11_inclusive", ";#lambda_{1}^{1};Entries", {HistType::kTH1F, {{100, 0., 1.}}}},
      {"h_lambda21_inclusive", ";#lambda_{2}^{1};Entries", {HistType::kTH1F, {{100, 0., 1.}}}},

      // Data histograms
      {"h_jet_pt_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_data", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_data", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_mass_data", ";m_{D_{S}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.7, 2.15}}}},
      {"h_ds_pt_data", ";#it{p}_{T,D_{S}} (GeV/#it{c});entries", {HistType::kTH1F, {{250, 0., 100.}}}},
      {"h_ds_eta_data", ";#eta_{D_{S}};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_ds_phi_data", ";#phi_{D_{S}};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

      {"h_ds_jet_projection_data", ";z^{D_{S},jet}_{||};entries", {HistType::kTH1F, {{200, 0., 1.2}}}},
      {"h_ds_jet_distance_data", ";#DeltaR_{D_{S},jet};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
      {"h_ds_jet_mass_data", ";m_{jet}^{ch} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 0., 25.}}}},
      {"h_ds_jet_lambda11_data", ";#lambda_{1}^{1};entries", {HistType::kTH1F, {{100, 0., 1.0}}}},
      {"h_ds_jet_lambda21_data", ";#lambda_{2}^{1};entries", {HistType::kTH1F, {{100, 0., 1.0}}}},

      {"hSparse_ds_data", ";m_{D_{S}};#it{p}_{T,D_{S}};#it{p}_{T,jet};z^{D_{S},jet}_{||};#DeltaR_{D_{S},jet}", {HistType::kTHnSparseF, {{60, 1.6, 2.3}, {60, 0., 80.}, {60, 0., 100.}, {20, 0., 1.2}, {20, 0., 1.0}}}},

      // MC general histograms
      {"McEffJet", "N_{jet};", {HistType::kTH1F, {{4, 0., 4.0}}}},
      {"McEffCol", "N_{collisions};", {HistType::kTH1F, {{4, 0., 4.0}}}},

      // MC detector-level histograms
      {"h_jet_pt_mcd", "detector-level jet pT;#it{p}_{T,jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_mcd", "detector-level jet #eta;#eta_{jet}^{det};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_mcd", "detector-level jet #phi;#phi_{jet}^{det};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_pt_mcd", ";#it{p}_{T,D_{S} jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}}},
      {"h_ds_eta_mcd", ";#eta_{D_{S} jet}^{det};entries", {HistType::kTH1F, {{60, -1., 1.}}}},
      {"h_ds_phi_mcd", ";#phi_{D_{S} jet}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}}},
      {"h_ds_mass_mcd", ";m_{D_{S}}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 1.7, 2.15}}}},

      // MCD - Sparse 1: mass, p_{T,Ds}, p_{T,jet}, z|| and prompt/non-prompt
      {"hSparse_ds_mcd1", ";m_{D_{S}}^{rec};#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||, det};Origin(D_{S})", {HistType::kTHnSparseF, {{60, 1.6, 2.3}, {60, 0., 80.}, {60, 0., 100.}, {20, 0., 1.2}, {2, -0.5, 1.5}}}},

      // MCD - Sparse 2: p_{T,Ds}, p_{T,jet}, and DeltaR
      {"hSparse_ds_mcd2", ";#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};#DeltaR_{D_{S},jet}^{det}", {HistType::kTHnSparseF, {{60, 0., 80.}, {60, 0., 100.}, {20, 0., 1.0}}}},
      // MCD - Sparse 3: p_{T,jet}, z|| and DeltaR
      {"hSparse_ds_mcd3", ";#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||, det};#DeltaR_{D_{S},jet}^{det}", {HistType::kTHnSparseF, {{60, 0., 100.}, {20, 0., 1.2}, {20, 0., 1.0}}}},

      // MC particle-level histograms
      {"h_jet_pt_mcp", "particle-level jet pT;#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta_mcp", "particle-level jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi_mcp", "particle-level jet #phi;#phi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},

      {"h_ds_pt_mcp", ";#it{p}_{T,D_{S} jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}}},
      {"h_ds_eta_mcp", ";#eta_{D_{S} jet}^{part};entries", {HistType::kTH1F, {{60, -1., 1.}}}},
      {"h_ds_phi_mcp", ";#phi_{D_{S} jet}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

      // MCP - Sparse: p_{T,Ds}, p_{T,jet}, z|| and DeltaR
      {"hSparse_ds_mcp", ";#it{p}_{T,D_{S}}^{part};#it{p}_{T,jet}^{part};z^{D_{S},jet}_{||, part};#DeltaR_{D_{S},jet}^{part};Origin(D_{S};Matching status)", {HistType::kTHnSparseF, {{60, 0., 80.}, {60, 0., 100.}, {20, 0., 1.2}, {20, 0., 1.0}, {2, -0.5, 1.5}, {2, -0.5, 1.5}}}},
    }};
  //========
  // INIT
  //========

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto hEvt = registry.get<TH1>(HIST("h_event_counter_data"));
    hEvt->GetXaxis()->SetBinLabel(1, "Input collisions");
    hEvt->GetXaxis()->SetBinLabel(2, "Event selection");
    hEvt->GetXaxis()->SetBinLabel(3, "|z| < 10 cm");

    // Labels
    auto mcCollisionCounter = registry.get<TH1>(HIST("McEffCol"));
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::All, "mccollisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::ZCut, "z_cut");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::Matched, "collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::MatchedSel8ZCut, "sel8");

    auto jetCounter = registry.get<TH1>(HIST("McEffJet"));
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetInMCCollision, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetInMCCollision, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetWithMatchedCandidate, "detector matched jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetWithMatchedCandidate, "particle matched jets");

    auto hSparse_ds_mcd1 = registry.get<THnSparse>(HIST("hSparse_ds_mcd1"));
    auto* axisOrigin = hSparse_ds_mcd1->GetAxis(4);
    axisOrigin->SetBinLabel(1, "Prompt");
    axisOrigin->SetBinLabel(2, "Non-prompt");

    auto hSparse_ds_mcp = registry.get<THnSparse>(HIST("hSparse_ds_mcp"));
    auto* axisOriginMCP = hSparse_ds_mcp->GetAxis(4);
    axisOriginMCP->SetBinLabel(1, "Prompt");
    axisOriginMCP->SetBinLabel(2, "Non-prompt");

    auto* axisMatching = hSparse_ds_mcp->GetAxis(5);
    axisMatching->SetBinLabel(1, "Unmatched");
    axisMatching->SetBinLabel(2, "Matched");
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
    const float jetRadius = jet.r() / 100.f;

    const float denom = std::pow(jet.pt(), kappa) * std::pow(jetRadius, alpha);
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

  // This function is shared by the pp and Pb–Pb (event-wise subtraction)
  template <typename JetTable>
  void fillInclusiveJetQA(aod::JetCollision const& collision,
                          aod::JetTracks const& tracks,
                          JetTable const& jets)
  {
    // Collision counter
    registry.fill(HIST("h_collisions"), getValFromBin(AllCollisions));

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), getValFromBin(Sel8ZCut));

    // Track QA
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }

    // Inclusive charged-jet QA
    registry.fill(HIST("h_njets_inclusive"), jets.size());

    for (auto const& jet : jets) {

      registry.fill(HIST("h_jet_pt_inclusive"), jet.pt());
      registry.fill(HIST("h_jet_eta_inclusive"), jet.eta());
      registry.fill(HIST("h_jet_phi_inclusive"), jet.phi());

      auto constituents = jet.template tracks_as<aod::JetTracks>();

      registry.fill(HIST("h_jet_nconst_inclusive"), constituents.size());

      registry.fill(HIST("h_jet_mass_inclusive"), computeJetMass(constituents));
      registry.fill(HIST("h_lambda11_inclusive"), computeLambda(jet, constituents, 1.f, 1.f));
      registry.fill(HIST("h_lambda21_inclusive"), computeLambda(jet, constituents, 2.f, 1.f));
    }
  }
  // Inclusive charged-jet QA for pp collisions
  void processCollisions(aod::JetCollision const& collision,
                         aod::JetTracks const& tracks,
                         FilteredChargedJets const& jets)
  {
    fillInclusiveJetQA(collision, tracks, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "collision QA", false);
  // Inclusive charged-jet QA for Pb–Pb collisions using
  // event-wise constituent-subtracted jets
  void processCollisionsEWS(aod::JetCollision const& collision,
                            aod::JetTracks const& tracks,
                            FilteredChargedJetsEWS const& jets)
  {
    fillInclusiveJetQA(collision, tracks, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processCollisionsEWS, "collision QA EWS", false);

  //=====================================================================================
  // DATA function
  //=====================================================================================

  // Common implementation of the Ds-in-jet analysis for data
  // Shared by the pp and Pb–Pb (event-wise subtraction) workflows
  template <typename JetTable>
  void analyseDataDsJet(aod::JetCollision const& collision,
                        JetTable const& jets)
  {
    registry.fill(HIST("h_event_counter_data"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter_data"), 2);

    if (std::abs(collision.posZ()) >= vertexZCut) {
      return;
    }
    registry.fill(HIST("h_event_counter_data"), 3);

    for (const auto& jet : jets) {

      registry.fill(HIST("h_jet_pt_data"), jet.pt());
      registry.fill(HIST("h_jet_eta_data"), jet.eta());
      registry.fill(HIST("h_jet_phi_data"), jet.phi());

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);
      const float lambda21 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);

      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Loop over Ds candidates
      for (const auto& dsCandidate : jet.template candidates_as<aod::CandidatesDsData>()) {

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

      if (!jet.template candidates_as<aod::CandidatesDsData>().empty()) {
        // Jet Mass
        registry.fill(HIST("h_ds_jet_mass_data"), mjet);

        // Jet angularity
        if (lambda11 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda11_data"), lambda11);
        }
        if (lambda21 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda21_data"), lambda21);
        }
      }
    }
  }

  //=====================================================================================
  // DATA process
  //=====================================================================================

  // Data analysis using standard charged jets (pp)
  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      FilteredDsDataJets const& jets,
                                      aod::CandidatesDsData const&,
                                      aod::JetTracks const&)
  {
    analyseDataDsJet(collision, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "Data charged jets", false);

  // Data analysis using event-wise constituent-subtracted charged jets (Pb–Pb)
  void processDataChargedSubstructureEWS(aod::JetCollision const& collision,
                                         FilteredDsDataJetsEWS const& jets,
                                         aod::CandidatesDsData const&,
                                         aod::JetTracks const&)
  {
    analyseDataDsJet(collision, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructureEWS, "Data charged jets EWS", false);

  //=====================================================================================
  //  MC function
  //=====================================================================================
  template <typename MCDJetsPerCollisionPreslice,
            typename MCPJetsPerCollisionPreslice,
            typename MCDJetTable,
            typename MCPJetTable,
            typename DsCandidatesMCD,
            typename DsCandidatesMCP>
  void analyseMonteCarloEfficiency(MCDJetsPerCollisionPreslice const& jetmcdpreslice,
                                   MCPJetsPerCollisionPreslice const& jetmcppreslice,
                                   aod::JetMcCollisions const& mccollisions,
                                   aod::JetCollisionsMCD const& collisions,
                                   MCDJetTable const& mcdjets,
                                   MCPJetTable const& mcpjets,
                                   DsCandidatesMCD const&,
                                   DsCandidatesMCP const&)
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
          auto mcdDscand = mcdjet.template candidates_first_as<DsCandidatesMCD>();

          // Check if it's prompt
          int origin = (mcdDscand.originMcRec() != RecoDecay::OriginType::Prompt) ? 1 : 0;

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

          // Detector-level Jet Histograms
          registry.fill(HIST("h_jet_pt_mcd"), mcdjet.pt());
          registry.fill(HIST("h_jet_eta_mcd"), mcdjet.eta());
          registry.fill(HIST("h_jet_phi_mcd"), mcdjet.phi());

          // Detector-level Ds Histgrams
          registry.fill(HIST("h_ds_pt_mcd"), mcdDscand.pt());
          registry.fill(HIST("h_ds_mass_mcd"), mcdDscand.m());
          registry.fill(HIST("h_ds_eta_mcd"), mcdDscand.eta());
          registry.fill(HIST("h_ds_phi_mcd"), mcdDscand.phi());

          // MCD THnSparse1: invariant mass, p{T,Ds}, pT, z, and origin (prompt/non-prompt)
          registry.fill(HIST("hSparse_ds_mcd1"),
                        mcdDscand.m(),
                        mcdDscand.pt(),
                        mcdjet.pt(),
                        mcd_zParallel,
                        origin);
          // MCD THnSparse2: invariant p{T,Ds}, pT and DeltaR
          registry.fill(HIST("hSparse_ds_mcd2"),
                        mcdDscand.pt(),
                        mcdjet.pt(),
                        mcd_deltaR);
          // MCD THnSparse3: invariant pT z and DeltaR
          registry.fill(HIST("hSparse_ds_mcd3"),
                        mcdjet.pt(),
                        mcd_zParallel,
                        mcd_deltaR);
        }
      }
      // Particle level
      const auto dsmcpJetsPerMCCollision = mcpjets.sliceBy(jetmcppreslice, mccollision.globalIndex());
      for (const auto& mcpjet : dsmcpJetsPerMCCollision) {

        registry.fill(HIST("McEffJet"), getValFromBin(BinMCJetCntr::ParticleLevelJetInMCCollision));

        // obtain leading HF particle in jet
        auto mcpDscand = mcpjet.template candidates_first_as<DsCandidatesMCP>();

        // Check if it's prompt 0 prompt, 1 non-prompt
        int originMCP = (mcpDscand.originMcGen() != RecoDecay::OriginType::Prompt) ? 1 : 0;

        // Matching status: 1 if the particle-level jet has a detector-level partner, 0 otherwise
        int isMatched = mcpjet.has_matchedJetCand() ? 1 : 0;

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("McEffJet"), getValFromBin(BinMCJetCntr::ParticleLevelJetWithMatchedCandidate));
        }

        TVector3 mcp_jetvector(mcpjet.px(), mcpjet.py(), mcpjet.pz());
        TVector3 mcp_candvector(mcpDscand.px(), mcpDscand.py(), mcpDscand.pz());

        float mcp_zParallel = (mcp_jetvector * mcp_candvector) / (mcp_jetvector * mcp_jetvector);
        // Axis distance Delta_R
        float mcp_deltaR = jetutilities::deltaR(mcpjet, mcpDscand);

        // Particle-level Jet Histograms
        registry.fill(HIST("h_jet_pt_mcp"), mcpjet.pt());
        registry.fill(HIST("h_jet_eta_mcp"), mcpjet.eta());
        registry.fill(HIST("h_jet_phi_mcp"), mcpjet.phi());

        // Particle-level Ds Histgrams
        registry.fill(HIST("h_ds_pt_mcp"), mcpDscand.pt());
        registry.fill(HIST("h_ds_eta_mcp"), mcpDscand.eta());
        registry.fill(HIST("h_ds_phi_mcp"), mcpDscand.phi());

        // Main THnSparse: invariant mass, pT, z, DeltaR , and origin (prompt/non-prompt)
        registry.fill(HIST("hSparse_ds_mcp"),
                      mcpDscand.pt(),
                      mcpjet.pt(),
                      mcp_zParallel,
                      mcp_deltaR,
                      originMCP,
                      isMatched);
      }
    }
  }
  //=====================================================================================
  // MC process
  //=====================================================================================

  // MC efficiency analysis using standard Ds-tagged jets (pp)
  void processMonteCarloEfficiencyDs(aod::JetMcCollisions const& mccollisions,
                                     aod::JetCollisionsMCD const& collisions,
                                     FilteredDsMCDJets const& mcdjets,
                                     FilteredDsMCPJets const& mcpjets,
                                     DsCandidatesMCD const& mcdDscand,
                                     DsCandidatesMCP const& mcpDscand)
  {
    analyseMonteCarloEfficiency<Preslice<DsMCDJets>,
                                Preslice<DsMCPJets>,
                                FilteredDsMCDJets,
                                FilteredDsMCPJets,
                                DsCandidatesMCD,
                                DsCandidatesMCP>(dsMCDJetsPerEXPCollisionPreslice,
                                                 dsMCPJetsPerMCCollisionPreslice,
                                                 mccollisions,
                                                 collisions,
                                                 mcdjets,
                                                 mcpjets,
                                                 mcdDscand,
                                                 mcpDscand);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processMonteCarloEfficiencyDs, "Non-matched and matched MC Ds and jets", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc)};
}
