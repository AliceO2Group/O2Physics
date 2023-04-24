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
///
/// \brief this is a starting point for the Strangeness tutorial
/// \author
/// \since 12/05/2023
/// \file strangeness_step3.cxx
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of V0s
// STEP 2
// Apply PID selections on V0 daughter tracks
// STEP 3
// Check the MC information of the V0s and verify with the PID information of daughter tracks

struct strangeness_tutorial {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurables parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};

  AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}; // Definition of axis

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                              {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}}},
                              {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}}},
                              {"hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}}},
                              {"hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}}},
                              {"hMassK0ShortTruePions", "hMassK0ShortTruePions", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hMassK0ShortMCTrue", "hMassK0ShortMCTrue", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hPtK0ShortTrue", "hPtK0ShortTrue", {HistType::kTH1F, {{ptAxis}}}}}};

  // Defining filters for events (event selection)
  // Processed events will be already fulfulling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daughters only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv&& nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv&& aod::v0data::dcaV0daughters < v0setting_dcav0dau;

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::McTrackLabels>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s,
               DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
               aod::McParticles const& mcParticles)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());

    for (auto& v0 : V0s) {

      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();

      registry.fill(HIST("hMassK0Short"), v0.mK0Short());

      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0setting_cospa)
        continue;
      if (v0.v0radius() < v0setting_radius)
        continue;
      if (TMath::Abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      registry.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
      registry.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (0.45 > v0.mK0Short() && v0.mK0Short() < 0.55) {
        registry.fill(HIST("hNSigmaPosPionFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
        registry.fill(HIST("hNSigmaNegPionFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
      }

      if (posDaughterTrack.has_mcParticle() && negDaughterTrack.has_mcParticle()) { // Checking that the daughter tracks come from particles and are not fake
        auto posParticle = posDaughterTrack.mcParticle();
        auto negParticle = negDaughterTrack.mcParticle();
        if (posParticle.pdgCode() == 211 && negParticle.pdgCode() == -211) { // Checking that the daughter tracks are true pions
          registry.fill(HIST("hMassK0ShortTruePions"), v0.mK0Short());
        }
      }

      // Checking that the V0 is a true K0s in the MC
      if (v0.has_mcParticle()) {
        auto v0mcparticle = v0.mcParticle();
        if (v0mcparticle.pdgCode() == 310) {
          registry.fill(HIST("hMassK0ShortMCTrue"), v0.mK0Short());
        }
      }
    }

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == 310) {
        registry.fill(HIST("hPtK0ShortTrue"), mcParticle.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)}; }
