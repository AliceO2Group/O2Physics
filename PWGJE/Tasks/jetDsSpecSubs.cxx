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

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "Common/Core/RecoDecay.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include <cmath>
#include <vector>
#include <string>

#include <TH1.h>
#include <TVector3.h>
#include <Math/Vector4D.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDsSpecSubs {

  //==================
  // Type definitions
  //==================

  using DsCandidatesData = aod::CandidatesDsData;
  using DsCandidatesMCD = aod::CandidatesDsMCD;
  using DsCandidatesMCP = aod::CandidatesDsMCP;

  using DsDataJets = soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents>;
  using DsMCDJets = soa::Join<aod::DsChargedMCDetectorLevelJets, aod::DsChargedMCDetectorLevelJetConstituents>;
  using DsMCPJets = soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents>;

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
  Filter jetCuts = aod::jet::pt > jetPtMin && aod::jet::r == nround(jetR.node() * 100.0f);
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  //=============
  // Histograms
  //=============

  HistogramRegistry registry{
    "registry",
    {

      {"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}},
      {"h_collision_counter", ";event counter;entries", {HistType::kTH1F, {{10, 0., 10.}}}},

      {"h_track_pt", ";#it{p}_{T,track};entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

      {"h_jet_counter", ";type;counts", {HistType::kTH1F, {{2, 0., 2.}}}},

      {"h_jet_pt", ";#it{p}_{T,jet};entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta", ";#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_jet_phi", ";#phi_{jet};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

      {"h_ds_mass", ";m_{D_{S}};entries", {HistType::kTH1F, {{400, 1.7, 2.2}}}},
      {"h_ds_pt", ";p_{T,D_{S}};entries", {HistType::kTH1F, {{200, 0., 100.}}}},
      {"h_ds_eta", ";#eta_{D_{S}};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"h_ds_phi", ";#phi_{D_{S}};entries", {HistType::kTH1F, {{100, -1., 7.}}}},

      {"h_ds_jet_projection", ";z_{||}^{D_{S},jet};entries", {HistType::kTH1F, {{200, 0., 2.}}}},
      {"h_ds_jet_distance", ";#DeltaR_{D_{S},jet};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
      {"h_ds_jet_lambda11", ";#lambda_{1}^{1};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
      {"h_ds_jet_lambda12", ";#lambda_{2}^{1};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
      {"h_ds_jet_mass", ";m_{jet};entries", {HistType::kTH1F, {{200, 0., 50.}}}},

      {"hSparse_ds", ";m_{D_{S}};p_{T,D_{S}};p_{T,jet};z_{||};#DeltaR",
       {HistType::kTHnSparseF,
        {
          {200, 1.7, 2.2},
          {200, 0., 100.},
          {200, 0., 100.},
          {200, 0., 2.},
          {200, 0., 1.}
        }}},

      {"h2_response_jet_pt", ";p_{T}^{det};p_{T}^{part}",
       {HistType::kTH2F,
        {
          {200, 0., 100.},
          {200, 0., 100.}
        }}},

      {"h2_response_lambda11", ";#lambda_{1}^{1,det};#lambda_{1}^{1,part}",
       {HistType::kTH2F,
        {
          {200, 0., 1.},
          {200, 0., 1.}
        }}}
    }
  };

  //========
  // INIT
  //========

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto h = registry.get<TH1>(HIST("h_jet_counter"));

    h->GetXaxis()->SetBinLabel(1, "All jets");
    h->GetXaxis()->SetBinLabel(2, "Ds-tagged jets");
  }

  //===============
  // Lambda compute
  //===============

  template<typename JET, typename TRACKS>
  float computeLambda(JET const& jet,
                      TRACKS const& tracks,
                      float alpha,
                      float kappa)
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

  template<typename TRACKS>
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

  void processCollisions(aod::JetCollision const& collision,
                         aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_collisions"), 1.5);

    for (auto const& track : tracks) {

      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }

  PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "collision QA", false);

  //===============
  // DATA analysis
  //===============

  template<typename JETS, typename CANDS>
  void analyseData(JETS const& jets,
                   CANDS const&)
  {
    for (const auto& jet : jets) {

      registry.fill(HIST("h_jet_counter"), 0.5);

      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);

      const float lambda12 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);

      bool hasDs = false;

      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (const auto& ds : jet.template candidates_as<CANDS>()) {

        hasDs = true;

        TVector3 dsVector(ds.px(), ds.py(), ds.pz());

        const float deltaR = jetutilities::deltaR(jet, ds);

        const float zParallel = (jetVector * dsVector) / (jetVector * jetVector);

        registry.fill(HIST("h_ds_mass"), ds.m());
        registry.fill(HIST("h_ds_pt"), ds.pt());
        registry.fill(HIST("h_ds_eta"), ds.eta());
        registry.fill(HIST("h_ds_phi"), ds.phi());

        registry.fill(HIST("h_ds_jet_projection"), zParallel);

        registry.fill(HIST("h_ds_jet_distance"), deltaR);

        registry.fill(HIST("hSparse_ds"),
                      ds.m(),
                      ds.pt(),
                      jet.pt(),
                      zParallel,
                      deltaR);
      }

      if (hasDs) {

        registry.fill(HIST("h_jet_counter"), 1.5);

        registry.fill(HIST("h_ds_jet_lambda11"), lambda11);
        registry.fill(HIST("h_ds_jet_lambda12"), lambda12);
        registry.fill(HIST("h_ds_jet_mass"), mjet);
      }
    }
  }

  //==============
  // MCD analysis
  //==============

  template<typename JETS>
  void analyseMCD(JETS const& jets)
  {
    for (const auto& jet : jets) {

      auto jetTracks =
        jet.template tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);

      const float lambda12 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);

      TVector3 jetVector(jet.px(), jet.py(),jet.pz());

      for (const auto& ds : jet.template candidates_as<DsCandidatesMCD>()) {

        TVector3 dsVector(ds.px(), ds.py(), ds.pz());

        const float deltaR = jetutilities::deltaR(jet, ds);

        const float zParallel = (jetVector * dsVector) / (jetVector * jetVector);

        registry.fill(HIST("hSparse_ds"),
                      ds.m(),
                      ds.pt(),
                      jet.pt(),
                      zParallel,
                      deltaR);

        registry.fill(HIST("h_ds_jet_lambda11"), lambda11);
        registry.fill(HIST("h_ds_jet_lambda12"), lambda12);
        registry.fill(HIST("h_ds_jet_mass"), mjet);
      }
    }
  }

  //==================
  // MC particle level
  //==================

  template<typename JETS>
  void analyseMCP(JETS const& jets)
  {
    for (const auto& jet : jets) {

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);

      const float lambda12 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);

      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (const auto& ds : jet.template candidates_as<DsCandidatesMCP>()) {

        TVector3 dsVector(ds.px(), ds.py(), ds.pz());

        const float deltaR = jetutilities::deltaR(jet, ds);

        const float zParallel = (jetVector * dsVector) / (jetVector * jetVector);

        registry.fill(HIST("hSparse_ds"),
                      ds.pt(),
                      jet.pt(),
                      zParallel,
                      deltaR);

        registry.fill(HIST("h_ds_jet_lambda11"), lambda11);
        registry.fill(HIST("h_ds_jet_lambda12"), lambda12);
        registry.fill(HIST("h_ds_jet_mass"), mjet);
      }
    }
  }

  //==============
  // DATA process
  //==============

  void processDataChargedSubstructure(
    aod::JetCollision const& collision,
    soa::Filtered<DsDataJets> const& jets,
    DsCandidatesData const& candidates,
    aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    analyseData(jets, candidates);
  }

  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "data charged jets", false);

  //==============
  // MCD process
  //==============

  void processMCDChargedSubstructure(
    aod::JetCollision const& collision,
    soa::Filtered<DsMCDJets> const& jets,
    aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    analyseMCD(jets);
  }

  PROCESS_SWITCH(JetDsSpecSubs, processMCDChargedSubstructure, "MC detector level", false);

  //==============
  // MCP process
  //==============

  void processMCPChargedSubstructure(
    aod::JetCollision const& collision,
    soa::Filtered<DsMCPJets> const& jets,
    aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    analyseMCP(jets);
  }

  PROCESS_SWITCH(JetDsSpecSubs, processMCPChargedSubstructure, "MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc)};
}