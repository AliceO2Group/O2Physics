// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetDrift.cxx
/// \brief Refactored jet drift correlation analysis task using O2Physics framework
/// \author Seokhwan / refactored by Antigravity

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
static constexpr int nPhiRegions = 8;
static const char* const phiLabel[8] = {"Phi0", "PhiPi4", "PhiPi2", "Phi3Pi4",
                                        "PhiPi", "Phi5Pi4", "Phi3Pi2", "Phi7Pi4"};

struct TrackInfo {
  float eta;
  float phi;
  float pt;
};

struct BufferedEvent {
  std::vector<TrackInfo> particles[2]; // Binned by associate pt
  double psiEP;
};

static inline float DeltaPhi_mpi2_3pi2(float phi_trig, float phi_assoc)
{
  float dphi = phi_trig - phi_assoc;
  while (dphi > TMath::Pi()) {
    dphi -= 2.0f * TMath::Pi();
  }
  while (dphi <= -TMath::Pi()) {
    dphi += 2.0f * TMath::Pi();
  }
  if (dphi < -TMath::Pi() / 2.0f) {
    dphi += 2.0f * TMath::Pi();
  }
  return dphi;
}

static inline int JetPhiRegion(float phi, float psiEP)
{
  float dphi = phi - psiEP;
  while (dphi < -TMath::Pi() / 6.0) {
    dphi += 2.0 * TMath::Pi();
  }
  while (dphi >= 11.0 * TMath::Pi() / 6.0) {
    dphi -= 2.0 * TMath::Pi();
  }

  if (dphi >= -TMath::Pi() / 6.0 && dphi < TMath::Pi() / 6.0) {
    return 0;
  }
  if (dphi >= TMath::Pi() / 6.0 && dphi < TMath::Pi() / 3.0) {
    return 1;
  }
  if (dphi >= TMath::Pi() / 3.0 && dphi < 2.0 * TMath::Pi() / 3.0) {
    return 2;
  }
  if (dphi >= 2.0 * TMath::Pi() / 3.0 && dphi < 5.0 * TMath::Pi() / 6.0) {
    return 3;
  }
  if (dphi >= 5.0 * TMath::Pi() / 6.0 && dphi < 7.0 * TMath::Pi() / 6.0) {
    return 4;
  }
  if (dphi >= 7.0 * TMath::Pi() / 6.0 && dphi < 4.0 * TMath::Pi() / 3.0) {
    return 5;
  }
  if (dphi >= 4.0 * TMath::Pi() / 3.0 && dphi < 5.0 * TMath::Pi() / 3.0) {
    return 6;
  }
  if (dphi >= 5.0 * TMath::Pi() / 3.0 && dphi < 11.0 * TMath::Pi() / 6.0) {
    return 7;
  }
  return -1;
}
} // namespace

struct jetDrift {
  EventPlaneHelper helperEP;
  HistogramRegistry registry{"registry", {// Region 0: Phi0
                                          {"PhiRegion_Phi0/hS_r0_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi0/hS_r0_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi0/hM_v2_r0_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi0/hM_v2_r0_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 1: PhiPi4
                                          {"PhiRegion_PhiPi4/hS_r1_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi4/hS_r1_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi4/hM_v2_r1_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi4/hM_v2_r1_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 2: PhiPi2
                                          {"PhiRegion_PhiPi2/hS_r2_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi2/hS_r2_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi2/hM_v2_r2_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi2/hM_v2_r2_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 3: Phi3Pi4
                                          {"PhiRegion_Phi3Pi4/hS_r3_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi4/hS_r3_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi4/hM_v2_r3_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi4/hM_v2_r3_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 4: PhiPi
                                          {"PhiRegion_PhiPi/hS_r4_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi/hS_r4_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi/hM_v2_r4_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_PhiPi/hM_v2_r4_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 5: Phi5Pi4
                                          {"PhiRegion_Phi5Pi4/hS_r5_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi5Pi4/hS_r5_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi5Pi4/hM_v2_r5_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi5Pi4/hM_v2_r5_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 6: Phi3Pi2
                                          {"PhiRegion_Phi3Pi2/hS_r6_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi2/hS_r6_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi2/hM_v2_r6_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi3Pi2/hM_v2_r6_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Region 7: Phi7Pi4
                                          {"PhiRegion_Phi7Pi4/hS_r7_a0", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi7Pi4/hS_r7_a1", "Signal;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi7Pi4/hM_v2_r7_a0", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"PhiRegion_Phi7Pi4/hM_v2_r7_a1", "Mixed Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Register global backgrounds
                                          {"InclusiveBackground/hB_global_a0", "Global Flat Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},
                                          {"InclusiveBackground/hB_global_a1", "Global Flat Background;#Delta#eta;#Delta#phi", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}}},

                                          // Register normalization tracking histograms
                                          {"hTrigW", "Trigger weight sum per region;region;weight", {HistType::kTH1F, {{8, -0.5, 7.5}}}},

                                          // Register QA histograms
                                          {"QA/h_collisions", "Centrality;centrality (%);entries", {HistType::kTH1F, {{10, 0.0, 100.0}}}},
                                          {"QA/h_z_vertex", "Z-vertex distribution;z_{vtx} (cm);entries", {HistType::kTH1F, {{200, -20.0, 20.0}}}},
                                          {"QA/h_multiplicity", "Track multiplicity;N_{tracks};entries", {HistType::kTH1F, {{500, 0.0, 2000.0}}}},
                                          {"QA/h_q2_mag", "Q_{2} vector magnitude;|Q_{2}|;entries", {HistType::kTH1F, {{150, 0.0, 3.0}}}},
                                          {"QA/h_v2", "Event-averaged v_{2};#LTv_{2}#GT;entries", {HistType::kTH1F, {{200, -0.5f, 0.5f}}}},
                                          {"QA/p_v2_pt", "v_{2}(p_{T});p_{T} (GeV/c);v_{2}", {HistType::kTProfile, {{100, 0.0f, 10.0f}}}},
                                          {"QA/h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{2000, 0.0, 200.0}}}},
                                          {"QA/h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{1000, -1.0, 1.0}}}},
                                          {"QA/h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{800, -1.0, 7.0}}}},
                                          {"QA/h_jet_pt", "raw jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{2000, 0.0, 200.0}}}},
                                          {"QA/h_jet_pt_corr", "corrected jet pT;#it{p}_{T,jet,corr} (GeV/#it{c});entries", {HistType::kTH1F, {{2000, -50.0, 150.0}}}},
                                          {"QA/h_jet_eta", "raw jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{1000, -1.0, 1.0}}}},
                                          {"QA/h_jet_phi", "raw jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{800, -1.0, 7.0}}}}}};

  // Event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0f, "Minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0f, "Maximum centrality"};

  // Jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 50.0f, "Minimum jet pT"};
  Configurable<float> jetR{"jetR", 0.4f, "Jet resolution parameter"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.6f, "Maximum jet pseudorapidity (e.g. 1.0 - R)"};
  Configurable<float> areaMin{"areaMin", 0.35f, "Minimum jet area"};

  // Track level configurables
  Configurable<std::vector<float>> assocPtMin{"assocPtMin", {0.2f, 2.0f}, "Associate track pt minimum limits"};
  Configurable<std::vector<float>> assocPtMax{"assocPtMax", {2.0f, 4.0f}, "Associate track pt maximum limits"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "Choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "Set track selections"};

  // Mixing level configurables
  Configurable<int> numberEventsMixed{"numberEventsMixed", 100, "Number of events to mix"};
  Configurable<int> maxBufferSize{"maxBufferSize", 1000, "Maximum buffer size for mixed events"};

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  TRandom3 randGen{0};
  std::deque<BufferedEvent> eventBuffer;

  void fillSignal(int r, int a, double dEta, double dPhi)
  {
    if (a == 0) {
      switch (r) {
        case 0:
          registry.fill(HIST("PhiRegion_Phi0/hS_r0_a0"), dEta, dPhi);
          break;
        case 1:
          registry.fill(HIST("PhiRegion_PhiPi4/hS_r1_a0"), dEta, dPhi);
          break;
        case 2:
          registry.fill(HIST("PhiRegion_PhiPi2/hS_r2_a0"), dEta, dPhi);
          break;
        case 3:
          registry.fill(HIST("PhiRegion_Phi3Pi4/hS_r3_a0"), dEta, dPhi);
          break;
        case 4:
          registry.fill(HIST("PhiRegion_PhiPi/hS_r4_a0"), dEta, dPhi);
          break;
        case 5:
          registry.fill(HIST("PhiRegion_Phi5Pi4/hS_r5_a0"), dEta, dPhi);
          break;
        case 6:
          registry.fill(HIST("PhiRegion_Phi3Pi2/hS_r6_a0"), dEta, dPhi);
          break;
        case 7:
          registry.fill(HIST("PhiRegion_Phi7Pi4/hS_r7_a0"), dEta, dPhi);
          break;
      }
    } else {
      switch (r) {
        case 0:
          registry.fill(HIST("PhiRegion_Phi0/hS_r0_a1"), dEta, dPhi);
          break;
        case 1:
          registry.fill(HIST("PhiRegion_PhiPi4/hS_r1_a1"), dEta, dPhi);
          break;
        case 2:
          registry.fill(HIST("PhiRegion_PhiPi2/hS_r2_a1"), dEta, dPhi);
          break;
        case 3:
          registry.fill(HIST("PhiRegion_Phi3Pi4/hS_r3_a1"), dEta, dPhi);
          break;
        case 4:
          registry.fill(HIST("PhiRegion_PhiPi/hS_r4_a1"), dEta, dPhi);
          break;
        case 5:
          registry.fill(HIST("PhiRegion_Phi5Pi4/hS_r5_a1"), dEta, dPhi);
          break;
        case 6:
          registry.fill(HIST("PhiRegion_Phi3Pi2/hS_r6_a1"), dEta, dPhi);
          break;
        case 7:
          registry.fill(HIST("PhiRegion_Phi7Pi4/hS_r7_a1"), dEta, dPhi);
          break;
      }
    }
  }

  void fillMixed(int r, int a, double dEta, double dPhi, double weight)
  {
    if (a == 0) {
      switch (r) {
        case 0:
          registry.fill(HIST("PhiRegion_Phi0/hM_v2_r0_a0"), dEta, dPhi, weight);
          break;
        case 1:
          registry.fill(HIST("PhiRegion_PhiPi4/hM_v2_r1_a0"), dEta, dPhi, weight);
          break;
        case 2:
          registry.fill(HIST("PhiRegion_PhiPi2/hM_v2_r2_a0"), dEta, dPhi, weight);
          break;
        case 3:
          registry.fill(HIST("PhiRegion_Phi3Pi4/hM_v2_r3_a0"), dEta, dPhi, weight);
          break;
        case 4:
          registry.fill(HIST("PhiRegion_PhiPi/hM_v2_r4_a0"), dEta, dPhi, weight);
          break;
        case 5:
          registry.fill(HIST("PhiRegion_Phi5Pi4/hM_v2_r5_a0"), dEta, dPhi, weight);
          break;
        case 6:
          registry.fill(HIST("PhiRegion_Phi3Pi2/hM_v2_r6_a0"), dEta, dPhi, weight);
          break;
        case 7:
          registry.fill(HIST("PhiRegion_Phi7Pi4/hM_v2_r7_a0"), dEta, dPhi, weight);
          break;
      }
    } else {
      switch (r) {
        case 0:
          registry.fill(HIST("PhiRegion_Phi0/hM_v2_r0_a1"), dEta, dPhi, weight);
          break;
        case 1:
          registry.fill(HIST("PhiRegion_PhiPi4/hM_v2_r1_a1"), dEta, dPhi, weight);
          break;
        case 2:
          registry.fill(HIST("PhiRegion_PhiPi2/hM_v2_r2_a1"), dEta, dPhi, weight);
          break;
        case 3:
          registry.fill(HIST("PhiRegion_Phi3Pi4/hM_v2_r3_a1"), dEta, dPhi, weight);
          break;
        case 4:
          registry.fill(HIST("PhiRegion_PhiPi/hM_v2_r4_a1"), dEta, dPhi, weight);
          break;
        case 5:
          registry.fill(HIST("PhiRegion_Phi5Pi4/hM_v2_r5_a1"), dEta, dPhi, weight);
          break;
        case 6:
          registry.fill(HIST("PhiRegion_Phi3Pi2/hM_v2_r6_a1"), dEta, dPhi, weight);
          break;
        case 7:
          registry.fill(HIST("PhiRegion_Phi7Pi4/hM_v2_r7_a1"), dEta, dPhi, weight);
          break;
      }
    }
  }

  void fillBGlobal(int a, double dEta, double dPhiFlat, double weight)
  {
    if (a == 0) {
      registry.fill(HIST("InclusiveBackground/hB_global_a0"), dEta, dPhiFlat, weight);
    } else {
      registry.fill(HIST("InclusiveBackground/hB_global_a1"), dEta, dPhiFlat, weight);
    }
  }

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  void processData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::QvectorFT0CVecs>::iterator const& collision,
                   soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets,
                   aod::JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(*collision, eventSelectionBits, false, false)) {
      return;
    }

    // 1. Z-vertex distribution (Cut 통과 전 원본 분포 확인용)
    registry.fill(HIST("QA/h_z_vertex"), collision.posZ());

    if (std::abs(collision.posZ()) > vertexZCut) {
      return;
    }

    float cent = collision.centFT0M();
    registry.fill(HIST("QA/h_collisions"), cent);
    if (cent < centralityMin || cent >= centralityMax) {
      return;
    }

    // 2. Multiplicity 계산 및 Histogram 채우기
    int nSelectedTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nSelectedTracks++;
      }
    }
    registry.fill(HIST("QA/h_multiplicity"), static_cast<float>(nSelectedTracks));

    // 3. Event-wise v2 (Q2 magnitude) 계산 및 Histogram 채우기
    double q2Re = collision.qvecFT0CReVec()[0];
    double q2Im = collision.qvecFT0CImVec()[0];
    double q2Mag = std::sqrt(q2Re * q2Re + q2Im * q2Im);
    registry.fill(HIST("QA/h_q2_mag"), q2Mag);

    // Fill track QA
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("QA/h_track_pt"), track.pt());
      registry.fill(HIST("QA/h_track_eta"), track.eta());
      registry.fill(HIST("QA/h_track_phi"), track.phi());
    }

    // 1. Reconstruct Event Plane (psiReco) using FT0C Q-vectors
    double psiReco = helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], 2);

    // Track v2 QA calculation relative to Event Plane
    double sumV2 = 0.0;
    int nV2Tracks = 0;
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      if (std::abs(track.eta()) < 0.8f && track.pt() >= 0.2f) {
        double v2 = std::cos(2.0 * (track.phi() - psiReco));
        registry.fill(HIST("QA/p_v2_pt"), track.pt(), v2);
        sumV2 += v2;
        nV2Tracks++;
      }
    }
    if (nV2Tracks > 0) {
      registry.fill(HIST("QA/h_v2"), static_cast<float>(sumV2 / nV2Tracks));
    }

    // 2. Select and bin associate tracks
    std::vector<TrackInfo> curTracks[2];
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      float pt = track.pt();
      float eta = track.eta();
      float phi = track.phi();
      if (std::abs(eta) > 1.0f) {
        continue;
      }
      for (int a = 0; a < 2; ++a) {
        if (pt >= assocPtMin->at(a) && pt < assocPtMax->at(a)) {
          curTracks[a].push_back({eta, phi, pt});
        }
      }
    }

    // 3. Process jets and fill same event/mixed event correlations
    for (auto const& jet : jets) {
      registry.fill(HIST("QA/h_jet_pt"), jet.pt());
      registry.fill(HIST("QA/h_jet_eta"), jet.eta());
      registry.fill(HIST("QA/h_jet_phi"), jet.phi());

      double ptCorr = jet.pt() - collision.rho() * jet.area();
      registry.fill(HIST("QA/h_jet_pt_corr"), ptCorr);

      if (ptCorr <= jetPtMin || std::abs(jet.eta()) >= jetEtaMax || jet.area() <= areaMin) {
        continue;
      }

      // Find the highest pt constituent track (trigger track)
      float highestPt = -1.0;
      float triggerEta = jet.eta();
      float triggerPhi = jet.phi();
      for (auto& constituent : jet.tracks_as<aod::JetTracks>()) {
        if (constituent.pt() > highestPt) {
          highestPt = constituent.pt();
          triggerEta = constituent.eta();
          triggerPhi = constituent.phi();
        }
      }

      int r = JetPhiRegion(jet.phi(), psiReco);
      if (r >= 0) {
        registry.fill(HIST("hTrigW"), static_cast<float>(r));

        // Same Event Correlation
        for (int a = 0; a < 2; ++a) {
          for (const auto& p : curTracks[a]) {
            // Exclude the trigger track
            if (std::abs(p.pt - highestPt) < 1e-5 &&
                std::abs(p.eta - triggerEta) < 1e-5 &&
                std::abs(p.phi - triggerPhi) < 1e-5) {
              continue;
            }
            double dEta = triggerEta - p.eta;
            double dPhi = DeltaPhi_mpi2_3pi2(triggerPhi, p.phi);
            fillSignal(r, a, dEta, dPhi);
          }
        }

        // Mixed Event Correlation
        if (static_cast<int>(eventBuffer.size()) >= numberEventsMixed) {
          for (int a = 0; a < 2; ++a) {
            int nMix = curTracks[a].size() * 100;
            for (int im = 0; im < nMix; ++im) {
              const auto& bevt = eventBuffer[randGen.Integer(eventBuffer.size())];
              if (bevt.particles[a].empty()) {
                continue;
              }
              const auto& p = bevt.particles[a][randGen.Integer(bevt.particles[a].size())];

              // Flow Aligned Mixed Background
              double corrected_phi = p.phi - bevt.psiEP + psiReco;
              double dEta = triggerEta - p.eta;
              double dPhi = DeltaPhi_mpi2_3pi2(triggerPhi, corrected_phi);
              fillMixed(r, a, dEta, dPhi, 0.01);

              // Inclusive (Flat) Mixed Background
              double dPhiFlat = DeltaPhi_mpi2_3pi2(triggerPhi, p.phi);
              fillBGlobal(a, dEta, dPhiFlat, 0.01);
            }
          }
        }
      }
    }

    // 4. Update Event Buffer
    if (static_cast<int>(eventBuffer.size()) >= maxBufferSize) {
      eventBuffer.pop_front();
    }
    BufferedEvent bevt;
    bevt.particles[0] = curTracks[0];
    bevt.particles[1] = curTracks[1];
    bevt.psiEP = psiReco;
    eventBuffer.push_back(bevt);
  }

  PROCESS_SWITCH(jetDrift, processData, "Jet drift correlation task", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jetDrift>(cfgc, TaskName{"jet-drift"})};
}
