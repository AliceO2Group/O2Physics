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
/// \file strangeness_step4.cxx
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
// STEP Cascades
// Apply selections on topological variables of Cascades

struct strangeness_tutorial {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurables parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascradius"};
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}; // Definition of axis

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {
                               {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                               {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                               {"hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                               {"hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}}},
                               {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}}},
                               {"hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}}},
                               {"hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}}},
                               {"hMassXi", "hMassXi", {HistType::kTH1F, {{200, 1.28f, 1.36f}}}},
                               {"hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {{200, 1.28f, 1.36f}}}},
                               {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}}},
                               {"hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}}},
                             }};

  // Defining filters for events (event selection)
  // Processed events will be already fulfulling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daughters only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv&& nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv&& aod::v0data::dcaV0daughters < v0setting_dcav0dau;

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<aod::V0Datas> const& V0s,
               aod::CascDataExt const& Cascades, aod::V0sLinked const&,
               DaughterTracks const&)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());

    for (const auto& v0 : V0s) {

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
    }

    for (const auto& casc : Cascades) {

      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0Casc = v0index.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachDaughterTrackCasc = casc.bachelor_as<DaughterTracks>();
      auto posDaughterTrackCasc = v0Casc.posTrack_as<DaughterTracks>();
      auto negDaughterTrackCasc = v0Casc.negTrack_as<DaughterTracks>();

      registry.fill(HIST("hMassXi"), casc.mXi());

      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_cospa)
        continue;
      if (casc.cascradius() < cascadesetting_cascradius)
        continue;
      if (casc.sign() < 0) {
        if (TMath::Abs(posDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(negDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      } else {
        if (TMath::Abs(negDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(posDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      }
      if (TMath::Abs(bachDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      registry.fill(HIST("hMassXiSelected"), casc.mXi());
      registry.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
      registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)}; }
