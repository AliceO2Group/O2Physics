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

// Jet substructure and spectrum task for D_s mesons
//
// This task is used to reconstruct and analyse jets containing charged D_s
// mesons
//
/// \author Monalisa Melo <monalisa.melo@cern.ch>, Universidade de São Paulo
//

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "TVector3.h"

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jet_distance
{
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);

// extra
DECLARE_SOA_COLUMN(JetMass, jetMass, float);
DECLARE_SOA_COLUMN(JetGirth, jetGirth, float);
DECLARE_SOA_COLUMN(JetThrust, jetThrust, float);     // lambda_2^1
DECLARE_SOA_COLUMN(JetLambda11, jetLambda11, float); // lambda_1^1
} // namespace jet_distance

DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2,
                  jet_distance::JetMass,
                  jet_distance::JetGirth,
                  jet_distance::JetThrust,
                  jet_distance::JetLambda11);
} // namespace o2::aod

struct JetDsSpecSubs {
  HistogramRegistry registry{
    "registry",
    {
      {"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}},
      {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
      {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
      {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
      {"h_collision_counter", "# of collisions;", {HistType::kTH1F, {{200, 0., 200.}}}},
      {"h_jet_counter", ";type;counts", {HistType::kTH1F, {{2, 0., 2.}}}},
      {"h_ds_jet_projection", ";z^{D_{S},jet}_{||};dN/dz^{D_{S},jet}_{||}", {HistType::kTH1F, {{1000, 0., 2.}}}},
      {"h_ds_jet_distance_vs_projection", ";#DeltaR_{D_{S},jet};z^{D_{S},jet}_{||}", {HistType::kTH2F, {{1000, 0., 1.}, {1000, 0., 2.}}}},
      {"h_ds_jet_distance", ";#DeltaR_{D_{S},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 1.}}}},
      {"h_ds_jet_pt", ";p_{T,D_{S} jet};dN/dp_{T,D_{S} jet}", {HistType::kTH1F, {{1000, 0., 100.}}}},
      {"h_ds_jet_eta", ";#eta_{D_{S} jet};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_jet_phi", ";#phi_{D_{S} jet};entries", {HistType::kTH1F, {{250, -1., 7.}}}},
      {"hSparse_ds", ";m_{D_{S}};p_{T,D_{S}};z^{D_{S},jet}_{||};#DeltaR_{D_{S},jet}", {HistType::kTHnSparseF, {{60, 1.7, 2.1}, {60, 0., 100.}, {60, 0., 2.}, {60, 0., 1.0}}}},
      {"h_ds_mass", ";m_{D_{S}} (GeV/c^{2});entries", {HistType::kTH1F, {{1000, 0., 6.}}}},
      {"h_ds_eta", ";#eta_{D_{S}};entries", {HistType::kTH1F, {{250, -1., 1.}}}},
      {"h_ds_phi", ";#phi_{D_{S}};entries", {HistType::kTH1F, {{250, -1., 7.}}}},
      {"h_ds_jet_mass", ";m_{jet}^{ch} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
      {"h_ds_jet_lambda11", ";#lambda_{1}^{1};entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"h_ds_jet_lambda12", ";#lambda_{2}^{1} (thrust);entries", {HistType::kTH1F, {{200, 0., 1.0}}}},
      {"h_ds_jet_girth", ";g (#equiv #lambda_{1}^{1}R);entries", {HistType::kTH1F, {{200, 0., 0.5}}}},
      {"h2_dsjet_pt_lambda11", ";#it{p}_{T,jet} (GeV/#it{c});#lambda_{1}^{1}", {HistType::kTH2F, {{100, 0., 100.}, {200, 0., 1.0}}}},
      {"h2_dsjet_pt_lambda12", ";#it{p}_{T,jet} (GeV/#it{c});#lambda_{2}^{1}", {HistType::kTH2F, {{100, 0., 100.}, {200, 0., 1.0}}}},
      {"h2_dsjet_pt_mass", ";#it{p}_{T,jet} (GeV/#it{c});m_{jet}^{ch} (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, 0., 100.}, {200, 0., 50.0}}}},
      {"h2_dsjet_pt_girth", ";#it{p}_{T,jet} (GeV/#it{c});g", {HistType::kTH2F, {{100, 0., 100.}, {200, 0., 0.5}}}},
    }};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  Produces<aod::JetDistanceTable> distJetTable;

  template <typename JET, typename TRACKS>
  float computeLambda(JET const& jet, TRACKS const& tracks, float a, float k)
  {
    if (jet.pt() <= 0.f) {
      return -1.f;
    }
    float sum = 0.f;
    for (auto const& trk : tracks) {
      const float dr = jetutilities::deltaR(jet, trk);
      sum += std::pow(trk.pt(), k) * std::pow(dr, a);
    }
    const float R = jet.r() / 100.f;
    const float denom = std::pow(jet.pt(), k) * std::pow(R, a);
    if (denom <= 0.f) {
      return -1.f;
    }
    return sum / denom;
  }

  template <typename TRACKS>
  float computeJetMassFromTracksMass(TRACKS const& tracks)
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

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto h = registry.get<TH1>(HIST("h_jet_counter"));
    h->GetXaxis()->SetBinLabel(1, "All jets");
    h->GetXaxis()->SetBinLabel(2, "Ds-tagged jets");
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  void processCollisions(aod::JetCollision const& collision, aod::JetTracks const& tracks)
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
  PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "process JE collisions", false);

  void processDataCharged(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          soa::Filtered<aod::ChargedJets> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());
    }
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataCharged, "charged jets in data", false);

  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents> const& jets,
                                      aod::CandidatesDsData const&,
                                      aod::JetTracks const&)
  {
    registry.fill(HIST("h_collision_counter"), 2.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
        !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), 3.0);

    for (const auto& jet : jets) {

      registry.fill(HIST("h_jet_counter"), 0.5);

      bool hasDs = false;

      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Compute jet-level quantities once (independent of Ds candidates)
      auto jetTracks = jet.tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);
      const float lambda12 = computeLambda(jet, jetTracks, 2.f, 1.f);
      const float mjet = computeJetMassFromTracksMass(jetTracks);

      const float R = jet.r() / 100.f;
      const float girth = (lambda11 >= 0.f) ? (lambda11 * R) : -1.f;

      // Loop over Ds candidates (particle level)
      for (const auto& dsCandidate : jet.candidates_as<aod::CandidatesDsData>()) {

        hasDs = true;

        TVector3 dsVector(dsCandidate.px(), dsCandidate.py(), dsCandidate.pz());

        // zParallel defined as longitudinal momentum fraction along the jet axis
        const double zParallel = (jetVector * dsVector) / (jetVector * jetVector);
        const double axisDistance = jetutilities::deltaR(jet, dsCandidate);

        // --- Ds-level observables ---
        registry.fill(HIST("h_ds_jet_projection"), zParallel);
        registry.fill(HIST("h_ds_jet_distance_vs_projection"), axisDistance, zParallel);
        registry.fill(HIST("h_ds_jet_distance"), axisDistance);

        registry.fill(HIST("h_ds_mass"), dsCandidate.m());
        registry.fill(HIST("h_ds_eta"), dsCandidate.eta());
        registry.fill(HIST("h_ds_phi"), dsCandidate.phi());

        const float mass = dsCandidate.m();
        const float pt = dsCandidate.pt();
        const float z = zParallel;
        const float dR = axisDistance;

        // Main THnSparse: invariant mass, pT, z, and ΔR
        registry.fill(HIST("hSparse_ds"), mass, pt, z, dR);

        // --- output table ---
        auto scores = dsCandidate.mlScores();
        constexpr int kScore0 = 0;
        constexpr int kScore1 = 1;
        constexpr int kScore2 = 2;

        const float s0 = (scores.size() > kScore0) ? scores[kScore0] : -999.f;
        const float s1 = (scores.size() > kScore1) ? scores[kScore1] : -999.f;
        const float s2 = (scores.size() > kScore2) ? scores[kScore2] : -999.f;

        distJetTable(static_cast<float>(axisDistance),
                     jet.pt(), jet.eta(), jet.phi(),
                     static_cast<int>(jetTracks.size()),
                     dsCandidate.pt(), dsCandidate.eta(), dsCandidate.phi(),
                     dsCandidate.m(), dsCandidate.y(),
                     s0, s1, s2,
                     mjet, girth, lambda12, lambda11);
      }

      // Jet-level quantities (filled once per jet containing at least one Ds)
      if (hasDs) {

        registry.fill(HIST("h_jet_counter"), 1.5);

        // Jet properties
        registry.fill(HIST("h_ds_jet_pt"), jet.pt());
        registry.fill(HIST("h_ds_jet_eta"), jet.eta());
        registry.fill(HIST("h_ds_jet_phi"), jet.phi());

        // Jet substructure observables
        if (lambda11 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda11"), lambda11);
          registry.fill(HIST("h2_dsjet_pt_lambda11"), jet.pt(), lambda11);
        }

        if (lambda12 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda12"), lambda12);
          registry.fill(HIST("h2_dsjet_pt_lambda12"), jet.pt(), lambda12);
        }

        registry.fill(HIST("h_ds_jet_mass"), mjet);
        registry.fill(HIST("h2_dsjet_pt_mass"), jet.pt(), mjet);

        if (girth >= 0.f) {
          registry.fill(HIST("h_ds_jet_girth"), girth);
          registry.fill(HIST("h2_dsjet_pt_girth"), jet.pt(), girth);
        }
      }
    }
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "charged HF jet substructure", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc, TaskName{"jet-ds-spectrum-subs"})};
}
