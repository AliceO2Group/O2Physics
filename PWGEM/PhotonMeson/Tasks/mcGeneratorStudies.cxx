// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file mcGeneratorStudies.cxx
/// \brief Task that produces the generated pT spectrum of a given particle for MC studies based on on the fly MC simulations
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt

#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TString.h>

#include <cstdlib>

using namespace o2;
using namespace o2::framework;

struct MCGeneratorStudies {
  HistogramRegistry mHistManager{"MCGeneratorStudyHistograms"};
  Configurable<int> cfgSelectedParticleCode{"cfgSelectedParticleCode", 111, "PDG code of the particle to be investigated"};
  Configurable<float> cfgMaxZVertex{"cfgMaxZVertex", 10, "Maximum absolute z-vertex distance (cm)"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.8, "Maximum absolute rapidity of selected generated particles"};
  expressions::Filter zVertexFilter = aod::mccollision::posZ < cfgMaxZVertex && aod::mccollision::posZ > -cfgMaxZVertex;

  void init(InitContext const&)
  {
    mHistManager.add("Yield", Form("Yield of primary selected particles in all MC collisions (PDG code %d);#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}}", cfgSelectedParticleCode.value), HistType::kTH1D, {{250, 0, 25}});
    mHistManager.add("FullYield", Form("Yield of selected particles in all MC collisions (PDG code %d);#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}}", cfgSelectedParticleCode.value), HistType::kTH1D, {{250, 0, 25}});
    mHistManager.add("MCGenXSec", "Generated cross section = GetBinContent(1)/GetBinContent(2)", HistType::kTH1D, {{2, 0, 2}});
    auto hX = mHistManager.get<TH1>(HIST("MCGenXSec"));
    hX->GetXaxis()->SetBinLabel(1, "#sigma_{gen} (b)");
    hX->GetXaxis()->SetBinLabel(2, "nStreams");
    hX->SetBinContent(2, 1); // Unity to enable normalization by number of streams
    hX->SetBinError(2, 0);
    mHistManager.add("RapidityCut", "#Delta y = GetBinContent(1)/GetBinContent(2)", HistType::kTH1D, {{2, 0, 2}});
    auto hY = mHistManager.get<TH1>(HIST("RapidityCut"));
    hY->GetXaxis()->SetBinLabel(1, "2#Delta y");
    hY->GetXaxis()->SetBinLabel(2, "nStreams");
    hY->SetBinContent(1, 2 * cfgRapidityCut);
    hY->SetBinContent(2, 1); // Unity to enable normalization by number of streams
    hY->SetBinError(2, 0);
    mHistManager.add("EventWeight", "MC collision weight;Event weight;#bf{#it{N}_{MC Collisions}}", HistType::kTH1D, {{1000, 0, 10}});
    mHistManager.add("MCGenXSecDistribution", "Generated cross section distribution;#sigma_{gen} (b);#bf{#it{N}_{MC Collisions}}", HistType::kTH1D, {{20000, 0, 2}});
    mHistManager.add("MCGenXSecErrDistribution", "Generated cross section error distribution;#sigma_{gen}^{err} (b);#bf{#it{N}_{MC Collisions}}", HistType::kTH1D, {{20000, 0, 2}});
    mHistManager.add("EventWeightSum", "Sum of event weights", HistType::kTH1D, {{2, 0., 2.}}, true);
    auto hE = mHistManager.get<TH1>(HIST("EventWeightSum"));
    hE->GetXaxis()->SetBinLabel(1, "#Sigma w");
    hE->GetXaxis()->SetBinLabel(2, "N_{MC Collisions}");
    mHistManager.add("zVertex", "zVertex position;#bf{z (cm)};#bf{#it{N}_{MC Collisions}}", HistType::kTH1D, {{500, -25., 25.}}, true);
  }

  void process(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, aod::McParticles const& mcParticles, aod::HepMCXSections const& crossSections)
  {
    const double w = mcCollision.weight();
    mHistManager.fill(HIST("EventWeight"), w);
    mHistManager.fill(HIST("EventWeightSum"), 0.5, w); // Fill with event weight (used for normalization)
    mHistManager.fill(HIST("EventWeightSum"), 1.5, 1); // Fill with number of MC collisions
    mHistManager.fill(HIST("zVertex"), mcCollision.posZ(), w);
    if (crossSections.size() == 1) {
      auto const xs = crossSections.begin();
      const double sigmaGen = xs.xsectGen() * 1E-12;    // From pb to b
      const double sigmaGenErr = xs.xsectErr() * 1E-12; // From pb to b
      mHistManager.fill(HIST("MCGenXSecDistribution"), sigmaGen);
      mHistManager.fill(HIST("MCGenXSecErrDistribution"), sigmaGenErr);

      auto hX = mHistManager.get<TH1>(HIST("MCGenXSec"));
      hX->SetBinContent(1, sigmaGen); // Set bin content to latest generated cross section and uncertainty
      hX->SetBinError(1, sigmaGenErr);
    } else {
      LOG(warn) << "Multiple cross sections found: " << crossSections.size();
    }
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == cfgSelectedParticleCode && std::abs(mcParticle.y()) < cfgRapidityCut) {
        mHistManager.fill(HIST("FullYield"), mcParticle.pt(), w);
        bool isPrimary = mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator();
        bool isFromWD = (aod::pwgem::photonmeson::utils::mcutil::IsFromWD(mcCollision, mcParticle, mcParticles)) > 0;
        if (isPrimary && !isFromWD) {
          mHistManager.fill(HIST("Yield"), mcParticle.pt(), w);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context) { return WorkflowSpec{adaptAnalysisTask<MCGeneratorStudies>(context, TaskName{"mc-generator-studies"})}; }
