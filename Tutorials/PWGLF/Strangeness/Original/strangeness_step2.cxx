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
/// \brief Step2 of the Strangeness tutorial
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"

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

struct strangeness_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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

  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // K0s reconstruction
    // Mass
    rKzeroShort.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {K0ShortMassAxis}});

    // K0s topological/PID cuts
    rKzeroShort.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daughters only
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<aod::V0Datas> const& V0s,
               DaughterTracks const&)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

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
      rKzeroShort.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), v0.v0cosPA());

      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (0.45 < v0.mK0Short() && v0.mK0Short() < 0.55) {
        rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
        rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
