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
/// \brief Step4 of the V0 strangeness tutorial on derived data (K0S, Lambda, AntiLambda)
/// \author Sara Pucillo
/// based on the O2AT strangeness derived-data tutorial from:
/// \author Romain Schotter
/// and on the original codes from:
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histograms
// STEP 1
// Apply selections on topological variables of V0s
// STEP 2
// Apply TPC PID selections on V0 daughter tracks
// STEP 3
// Apply TOF PID selections on V0 daughter tracks (if info is available)
// STEP 4
// Check the MC information of the V0s

struct StrangenessDerivedV0Step4 {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0Short{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rGenParticles{"genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 0.5, "DCA V0 daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA positive to PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA negative to PV"};
  Configurable<float> v0setting_cospa{"v0setting_cospa", 0.995, "V0 CosPA"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "V0 radius"};
  Configurable<float> v0setting_rapidity{"v0setting_rapidity", 0.5, "V0 rapidity window"};

  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  // Configurable parameters for TOF PID selection
  Configurable<float> NSigmaTOFPion{"NSigmaTOFPion", 3, "NSigmaTOFPion"};
  Configurable<float> NSigmaTOFProton{"NSigmaTOFProton", 3, "NSigmaTOFProton"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {100, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {100, 1.08f, 1.15f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec nSigmaAxis = {100, -10.0f, 10.0f, "N#sigma_{TPC}"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // K0S/Lambda/AntiLambda reconstruction
    rK0Short.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
    rK0Short.add("hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {K0ShortMassAxis}});
    rK0Short.add("hMassK0ShortSelectedWithTOF", "hMassK0ShortSelectedWithTOF", {HistType::kTH1F, {K0ShortMassAxis}});
    rK0Short.add("hMassK0ShortTrueRec", "hMassK0ShortTrueRec", {HistType::kTH1F, {K0ShortMassAxis}});
    rK0Short.add("hPtK0ShortTrueRec", "hPtK0ShortTrueRec", {HistType::kTH1F, {ptAxis}});
    rK0Short.add("hMassK0ShortTrueRecWithTOF", "hMassK0ShortTrueRecWithTOF", {HistType::kTH1F, {K0ShortMassAxis}});
    rK0Short.add("hPtK0ShortTrueRecWithTOF", "hPtK0ShortTrueRecWithTOF", {HistType::kTH1F, {ptAxis}});

    rLambda.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelected", "hMassLambdaSelected", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambdaSelected", "hMassAntiLambdaSelected", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelectedWithTOF", "hMassLambdaSelectedWithTOF", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambdaSelectedWithTOF", "hMassAntiLambdaSelectedWithTOF", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaTrueRec", "hMassLambdaTrueRec", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hPtLambdaTrueRec", "hPtLambdaTrueRec", {HistType::kTH1F, {ptAxis}});
    rLambda.add("hMassLambdaTrueRecWithTOF", "hMassLambdaTrueRecWithTOF", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hPtLambdaTrueRecWithTOF", "hPtLambdaTrueRecWithTOF", {HistType::kTH1F, {ptAxis}});
    rLambda.add("hMassAntiLambdaTrueRec", "hMassAntiLambdaTrueRec", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hPtAntiLambdaTrueRec", "hPtAntiLambdaTrueRec", {HistType::kTH1F, {ptAxis}});
    rLambda.add("hMassAntiLambdaTrueRecWithTOF", "hMassAntiLambdaTrueRecWithTOF", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hPtAntiLambdaTrueRecWithTOF", "hPtAntiLambdaTrueRecWithTOF", {HistType::kTH1F, {ptAxis}});

    // V0 topological variables (after selections)
    rK0Short.add("hV0DCADaughters", "hV0DCADaughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rK0Short.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rK0Short.add("hV0Radius", "hV0Radius", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});

    // TPC PID of the K0S daughters, in the K0S peak region, before the PID selection
    rK0Short.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
    rK0Short.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {ptAxis, nSigmaAxis}});

    // Generated level histograms
    rEventSelection.add("hVertexZGen", "hVertexZGen", {HistType::kTH1F, {vertexZAxis}});
    rGenParticles.add("hPtK0ShortGen", "hPtK0ShortGen", {HistType::kTH1F, {ptAxis}});
    rGenParticles.add("hPtLambdaGen", "hPtLambdaGen", {HistType::kTH1F, {ptAxis}});
    rGenParticles.add("hPtAntiLambdaGen", "hPtAntiLambdaGen", {HistType::kTH1F, {ptAxis}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns
  Filter preFilterV0 = (aod::v0data::dcaV0daughters < v0setting_dcav0dau &&
                        nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv);

  // Defining the type of the daughter tracks
  using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

  void processRecMC(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision,
                    soa::Filtered<soa::Join<aod::V0Cores, aod::V0Extras, aod::V0TOFNSigmas, aod::V0CoreMCLabels>> const& V0s,
                    dauTracks const&,
                    aod::V0MCCores const& /*v0mccores*/)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // V0s
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrackExtra_as<dauTracks>();
      const auto& negDaughterTrack = v0.negTrackExtra_as<dauTracks>();

      rK0Short.fill(HIST("hMassK0Short"), v0.mK0Short());
      rLambda.fill(HIST("hMassLambda"), v0.mLambda());
      rLambda.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());

      // Keep only standard V0s (type 1): type 0 = not passing the standard V0 selections, 2-3 = photon-like (TPC-only)
      if (v0.v0Type() != 1) {
        continue;
      }
      // Topological selections not applied in the Filter
      // (in V0Cores the cosine of pointing angle is a regular column: it could also be moved into the Filter)
      if (v0.v0cosPA() < v0setting_cospa) {
        continue;
      }
      // Cut on dynamic columns
      if (v0.v0radius() < v0setting_radius) {
        continue;
      }
      // From here on, the same V0 is tested under the three mass hypotheses (K0S, Lambda, AntiLambda):
      // the hypothesis-dependent selections cannot use 'continue', they are stored in booleans

      // Kinematic selection
      bool isK0Short = std::abs(v0.yK0Short()) < v0setting_rapidity;
      bool isLambda = std::abs(v0.yLambda()) < v0setting_rapidity;
      bool isAntiLambda = std::abs(v0.yLambda()) < v0setting_rapidity;

      // TPC PID QA of the K0S daughters (in the K0S peak region, before the PID selection)
      if (isK0Short && std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < 0.01) {
        rK0Short.fill(HIST("hNSigmaPosPionFromK0s"), v0.positivept(), posDaughterTrack.tpcNSigmaPi());
        rK0Short.fill(HIST("hNSigmaNegPionFromK0s"), v0.negativept(), negDaughterTrack.tpcNSigmaPi());
      }

      // TPC PID selection
      // K0S -> pi+ pi-
      if (std::abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion || std::abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        isK0Short = false;
      }
      // Lambda -> p pi-
      if (std::abs(posDaughterTrack.tpcNSigmaPr()) > NSigmaTPCProton || std::abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        isLambda = false;
      }
      // AntiLambda -> pbar pi+
      if (std::abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion || std::abs(negDaughterTrack.tpcNSigmaPr()) > NSigmaTPCProton) {
        isAntiLambda = false;
      }

      // TOF PID check (applied only if the daughter has a TOF signal)
      bool k0ShortPassTOFSelection = true;
      bool lambdaPassTOFSelection = true;
      bool antiLambdaPassTOFSelection = true;
      if (v0.positiveHasTOF()) {
        if (std::abs(v0.tofNSigmaK0PiPlus()) > NSigmaTOFPion) {
          k0ShortPassTOFSelection &= false;
        }
        if (std::abs(v0.tofNSigmaLaPr()) > NSigmaTOFProton) {
          lambdaPassTOFSelection &= false;
        }
        if (std::abs(v0.tofNSigmaALaPi()) > NSigmaTOFPion) {
          antiLambdaPassTOFSelection &= false;
        }
      }
      if (v0.negativeHasTOF()) {
        if (std::abs(v0.tofNSigmaK0PiMinus()) > NSigmaTOFPion) {
          k0ShortPassTOFSelection &= false;
        }
        if (std::abs(v0.tofNSigmaLaPi()) > NSigmaTOFPion) {
          lambdaPassTOFSelection &= false;
        }
        if (std::abs(v0.tofNSigmaALaPr()) > NSigmaTOFProton) {
          antiLambdaPassTOFSelection &= false;
        }
      }
      // Equivalent short form: v0.tofK0ShortCompatibility(nSigma), v0.tofLambdaCompatibility(nSigma), ...

      // Fill histograms! (if possible)
      if (isK0Short) {
        rK0Short.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
        if (k0ShortPassTOFSelection) {
          rK0Short.fill(HIST("hMassK0ShortSelectedWithTOF"), v0.mK0Short());
        }

        rK0Short.fill(HIST("hV0DCADaughters"), v0.dcaV0daughters());
        rK0Short.fill(HIST("hV0CosPA"), v0.v0cosPA());
        rK0Short.fill(HIST("hV0Radius"), v0.v0radius());
      }
      if (isLambda) {
        rLambda.fill(HIST("hMassLambdaSelected"), v0.mLambda());
        if (lambdaPassTOFSelection) {
          rLambda.fill(HIST("hMassLambdaSelectedWithTOF"), v0.mLambda());
        }
      }
      if (isAntiLambda) {
        rLambda.fill(HIST("hMassAntiLambdaSelected"), v0.mAntiLambda());
        if (antiLambdaPassTOFSelection) {
          rLambda.fill(HIST("hMassAntiLambdaSelectedWithTOF"), v0.mAntiLambda());
        }
      }

      // MC truth info
      if (!v0.has_v0MCCore()) {
        continue;
      }
      auto v0mccore = v0.v0MCCore_as<aod::V0MCCores>();

      // Efficiency is computed for primary particles only (e.g. Lambda from Xi decays are excluded)
      if (!v0mccore.isPhysicalPrimary()) {
        continue;
      }

      // Checking that the V0 is a true K0S / Lambda / AntiLambda
      if (v0mccore.pdgCode() == PDG_t::kK0Short && isK0Short) {
        rK0Short.fill(HIST("hMassK0ShortTrueRec"), v0.mK0Short());
        rK0Short.fill(HIST("hPtK0ShortTrueRec"), v0.pt());
        if (k0ShortPassTOFSelection) {
          rK0Short.fill(HIST("hMassK0ShortTrueRecWithTOF"), v0.mK0Short());
          rK0Short.fill(HIST("hPtK0ShortTrueRecWithTOF"), v0.pt());
        }
      }
      if (v0mccore.pdgCode() == PDG_t::kLambda0 && isLambda) {
        rLambda.fill(HIST("hMassLambdaTrueRec"), v0.mLambda());
        rLambda.fill(HIST("hPtLambdaTrueRec"), v0.pt());
        if (lambdaPassTOFSelection) {
          rLambda.fill(HIST("hMassLambdaTrueRecWithTOF"), v0.mLambda());
          rLambda.fill(HIST("hPtLambdaTrueRecWithTOF"), v0.pt());
        }
      }
      if (v0mccore.pdgCode() == PDG_t::kLambda0Bar && isAntiLambda) {
        rLambda.fill(HIST("hMassAntiLambdaTrueRec"), v0.mAntiLambda());
        rLambda.fill(HIST("hPtAntiLambdaTrueRec"), v0.pt());
        if (antiLambdaPassTOFSelection) {
          rLambda.fill(HIST("hMassAntiLambdaTrueRecWithTOF"), v0.mAntiLambda());
          rLambda.fill(HIST("hPtAntiLambdaTrueRecWithTOF"), v0.pt());
        }
      }
    }
  }

  void processGenMC(soa::Filtered<aod::StraMCCollisions>::iterator const& mcCollision,
                    const soa::SmallGroups<soa::Join<aod::StraCollisions, o2::aod::StraEvSels, aod::StraCollLabels>>& collisions,
                    const soa::SmallGroups<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>& v0MC)
  {
    if (collisions.size() < 1) { // to process generated collisions that've been reconstructed at least once
      return;
    }
    rEventSelection.fill(HIST("hVertexZGen"), mcCollision.posZ());

    for (const auto& v0mc : v0MC) {
      if (!v0mc.isPhysicalPrimary()) {
        continue;
      }
      // rapidityMC(i): i = 0 -> K0S mass hypothesis, 1 -> Lambda, 2 -> AntiLambda
      if (v0mc.pdgCode() == PDG_t::kK0Short && std::abs(v0mc.rapidityMC(0)) < v0setting_rapidity) {
        rGenParticles.fill(HIST("hPtK0ShortGen"), v0mc.ptMC());
      }
      if (v0mc.pdgCode() == PDG_t::kLambda0 && std::abs(v0mc.rapidityMC(1)) < v0setting_rapidity) {
        rGenParticles.fill(HIST("hPtLambdaGen"), v0mc.ptMC());
      }
      if (v0mc.pdgCode() == PDG_t::kLambda0Bar && std::abs(v0mc.rapidityMC(2)) < v0setting_rapidity) {
        rGenParticles.fill(HIST("hPtAntiLambdaGen"), v0mc.ptMC());
      }
    }
  }

  PROCESS_SWITCH(StrangenessDerivedV0Step4, processRecMC, "Process Run 3 mc, reconstructed", true);
  PROCESS_SWITCH(StrangenessDerivedV0Step4, processGenMC, "Process Run 3 mc, generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessDerivedV0Step4>(cfgc, TaskName{"strangeness-derived-v0-tutorial"})};
}
