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
/// \brief Step4 of the Strangeness tutorial
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

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
// Check the MC information of the V0s
// STEP 4
// Apply selections on topological variables of Cascades

struct strangeness_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rXi{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rGenParticles{"genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};

  // Configurable parameters for cascade selection
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "DCA cascade daughters"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "DCA bachelor to PV"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "minimum V0 DCA to PV"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "v0 mass window"};

  Configurable<float> cascade_dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> cascade_dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};

  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec XiMassAxis = {200, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // K0s reconstruction
    // Mass
    rKzeroShort.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelectedTruePions", "hMassK0ShortSelectedTruePions", {HistType::kTH1F, {{200, 0.45f, 0.55f}}});
    rKzeroShort.add("hMassK0ShortTrueRec", "hMassK0ShortTrueRec", {HistType::kTH1F, {{200, 0.45f, 0.55f}}});
    // Pt
    rKzeroShort.add("hPtK0ShortSelected", "hPtK0ShortSelected", {HistType::kTH1F, {{ptAxis}}});
    rKzeroShort.add("hPtK0ShortTrueRec", "hPtK0ShortTrueRec", {HistType::kTH1F, {{ptAxis}}});

    // Xi reconstruction
    rXi.add("hMassXi", "hMassXi", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiTrueRec", "hMassXiTrueRec", {HistType::kTH1F, {XiMassAxis}});

    // K0s topological/PID cuts
    rKzeroShort.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});

    // Xi topological cuts
    rXi.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rXi.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});

    // Generated level histograms
    rEventSelection.add("hVertexZGen", "hVertexZGen", {HistType::kTH1F, {vertexZAxis}});
    rGenParticles.add("hPtK0ShortGen", "hPtK0ShortGen", {HistType::kTH1F, {{ptAxis}}});
    rGenParticles.add("hPtXiGen", "hPtXiGen", {HistType::kTH1F, {{ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Filters on Cascades
  Filter preFilterCascades = (nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascdata::dcaV0daughters < v0setting_dcav0dau &&
                              nabs(aod::cascdata::dcapostopv) > cascade_dcapostopv &&
                              nabs(aod::cascdata::dcanegtopv) > cascade_dcanegtopv &&
                              nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr, aod::McTrackLabels>;

  void processRecMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                    soa::Filtered<soa::Join<aod::CascDatas, aod::McCascLabels>> const& Cascades,
                    soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s,
                    DaughterTracks const&,
                    aod::McParticles const&)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // V0s
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();

      rKzeroShort.fill(HIST("hMassK0Short"), v0.mK0Short());

      // Cut on dynamic columns
      if (v0.v0cosPA() < v0setting_cospa)
        continue;
      if (v0.v0radius() < v0setting_radius)
        continue;

      if (std::abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (std::abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      rKzeroShort.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
      rKzeroShort.fill(HIST("hPtK0ShortSelected"), v0.pt());

      rKzeroShort.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), v0.v0cosPA());

      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (0.45 < v0.mK0Short() && v0.mK0Short() < 0.55) {
        rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
        rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
      }

      if (posDaughterTrack.has_mcParticle() && negDaughterTrack.has_mcParticle()) { // Checking that the daughter tracks come from particles and are not fake
        auto posParticle = posDaughterTrack.mcParticle();
        auto negParticle = negDaughterTrack.mcParticle();
        if (posParticle.pdgCode() == PDG_t::kPiPlus && negParticle.pdgCode() == PDG_t::kPiMinus) { // Checking that the daughter tracks are true pions
          rKzeroShort.fill(HIST("hMassK0ShortSelectedTruePions"), v0.mK0Short());
        }
      }

      // Checking that the V0 is a true K0s
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (v0mcParticle.pdgCode() == PDG_t::kK0Short) {
          rKzeroShort.fill(HIST("hMassK0ShortTrueRec"), v0.mK0Short());
          rKzeroShort.fill(HIST("hPtK0ShortTrueRec"), v0.pt()); // To mimic distribution after the signal extraction
        }
      }
    }

    // Cascades
    for (const auto& casc : Cascades) {
      const auto& bachDaughterTrackCasc = casc.bachelor_as<DaughterTracks>();
      const auto& posDaughterTrackCasc = casc.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrackCasc = casc.negTrack_as<DaughterTracks>();

      rXi.fill(HIST("hMassXi"), casc.mXi());

      // Cut on dynamic columns
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_cospa)
        continue;
      if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda) > cascadesetting_v0masswindow)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_mindcav0topv)
        continue;
      if (casc.cascradius() < cascadesetting_cascradius)
        continue;

      if (casc.sign() < 0) {
        if (std::abs(posDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (std::abs(negDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      } else {
        if (std::abs(negDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (std::abs(posDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      }
      if (std::abs(bachDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      rXi.fill(HIST("hMassXiSelected"), casc.mXi());
      rXi.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
      rXi.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));

      // Checking that the cascade is a true Xi
      if (casc.has_mcParticle()) {
        const auto cascmcParticle = casc.mcParticle();
        if (std::abs(cascmcParticle.pdgCode()) == PDG_t::kXiMinus) {
          rXi.fill(HIST("hMassXiTrueRec"), casc.mXi());
        }
      }
    }
  }

  void processGenMC(soa::Filtered<aod::McCollisions>::iterator const& mcCollision,
                    const soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>& collisions,
                    aod::McParticles const& mcParticles)
  {
    if (collisions.size() < 1) // to process generated collisions that've been reconstructed at least once
      return;
    rEventSelection.fill(HIST("hVertexZGen"), mcCollision.posZ());
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == PDG_t::kK0Short) {
        rGenParticles.fill(HIST("hPtK0ShortGen"), mcParticle.pt());
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus) {
        rGenParticles.fill(HIST("hPtXiGen"), mcParticle.pt());
      }
    }
  }

  PROCESS_SWITCH(strangeness_tutorial, processRecMC, "Process Run 3 mc, reconstructed", true);
  PROCESS_SWITCH(strangeness_tutorial, processGenMC, "Process Run 3 mc, generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
