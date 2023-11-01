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
/// \file strangeness_step1.cxx
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of V0s
struct strangeness_tutorial {

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurables parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{"registry",
                             {{"hVertexZ", "hVertexZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
                              {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {{200, 0.45f, 0.55f}}}},
                              {"hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}}},
                              {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}}}}};

  // Defining filters for events (event selection)
  // Processed events will be already fulfulling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daughters only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv&& nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv&& aod::v0data::dcaV0daughters < v0setting_dcav0dau;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<aod::V0Datas> const& V0s)
  {
    // Fill the event counter
    registry.fill(HIST("hVertexZ"), collision.posZ());

    for (auto& v0 : V0s) {
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());

      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0setting_cospa)
        continue;
      if (v0.v0radius() < v0setting_radius)
        continue;

      registry.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
      registry.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)}; }
