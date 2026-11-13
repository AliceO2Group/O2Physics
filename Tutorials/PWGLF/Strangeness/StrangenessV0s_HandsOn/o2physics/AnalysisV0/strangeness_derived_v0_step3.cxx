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
/// \brief Step3 of the V0 strangeness tutorial on derived data (K0S, Lambda, AntiLambda)
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

struct strangeness_derived_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0Short{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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

    rLambda.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelected", "hMassLambdaSelected", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambdaSelected", "hMassAntiLambdaSelected", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelectedWithTOF", "hMassLambdaSelectedWithTOF", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassAntiLambdaSelectedWithTOF", "hMassAntiLambdaSelectedWithTOF", {HistType::kTH1F, {LambdaMassAxis}});

    // V0 topological variables (after selections)
    rK0Short.add("hV0DCADaughters", "hV0DCADaughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rK0Short.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rK0Short.add("hV0Radius", "hV0Radius", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});

    // TPC PID of the K0S daughters, in the K0S peak region, before the PID selection
    rK0Short.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
    rK0Short.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {ptAxis, nSigmaAxis}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns
  Filter preFilterV0 = (aod::v0data::dcaV0daughters < v0setting_dcav0dau &&
                        nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv);

  // Defining the type of the daughter tracks
  using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

  void process(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision,
               soa::Filtered<soa::Join<aod::V0Cores, aod::V0Extras, aod::V0TOFNSigmas>> const& V0s,
               dauTracks const&)
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
      if (v0.v0Type() != 1)
        continue;

      // Topological selections not applied in the Filter
      // (in V0Cores the cosine of pointing angle is a regular column: it could also be moved into the Filter)
      if (v0.v0cosPA() < v0setting_cospa)
        continue;
      // Cut on dynamic columns
      if (v0.v0radius() < v0setting_radius)
        continue;

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
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_derived_tutorial>(cfgc)};
}
