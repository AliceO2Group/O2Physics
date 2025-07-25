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
/// \brief this is a starting point for the Resonances tutorial for MC
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>
/// \since 11/10/2024

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 4
// Producing histograms for Generated MCs
struct resonances_tutorial {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  /// Histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Binning of the centrality axis"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};

    // register histograms
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});

    // MC QA
    histos.add("hphipt", "pT distribution of True MC phi", kTH1F, {ptAxis});
    histos.add("phiGen", "pT distribution of True MC phi", kTH2F, {ptAxis, centAxis});

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in resonance tutorial step1:";
    histos.print();
  }

  // MC particle selection
  template <typename ParticleType>
  bool ptCut(const ParticleType resoParents)
  {
    // basic pt cuts
    if (std::abs(resoParents.pt()) < cMinPtcut)
      return false;

    return true;
  }

  // Fill histograms (main function)
  template <typename CollisionType, typename ParticleType>
  void fillHistograms(const CollisionType& collision, const ParticleType& resoParents)
  {
    auto multiplicity = collision.cent();
    for (auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (!ptCut(part))
        continue; // pt selection

      // QA plots
      histos.fill(HIST("hEta"), part.eta());

      if (abs(part.pdgCode()) != 333) // phi(0)
        continue;
      if (abs(part.y()) > 0.5) { // rapidity cut
        continue;
      }
      if (abs(part.daughterPDG1()) != 321 || abs(part.daughterPDG2()) != 321) { // At least one decay to Kaon
        continue;
      }
      histos.fill(HIST("phiGen"), part.pt(), multiplicity);
      histos.fill(HIST("hphipt"), part.pt());
    }
  }

  // Process the MC
  void process(aod::ResoCollision& collision, aod::ResoMCParents& resoParents)
  {
    auto multiplicity = collision.cent();

    // Fill the event counter
    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), multiplicity);

    fillHistograms(collision, resoParents); // Fill histograms, MC
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
