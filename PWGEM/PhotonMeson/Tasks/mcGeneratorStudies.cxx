// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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
/// \file mcGeneratorStudies.cxx
///
/// \brief Task that produces the generated pT spectrum of a given particle for MC studies based on on the fly MC simulations
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;

struct MCGeneratorStudies {
  HistogramRegistry mHistManager{"MCGeneratorStudyHistograms"};
  Configurable<int> cfgSelectedParticleCode{"cfgSelectedParticleCode", 111, "PDG code of the particle to be investigated"};
  Configurable<float> cfgMaxZVertex{"cfgMaxZVertex", 10, "Maximum absolute z-vertex distance (cm)"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.8, "Maximum absolute rapditity of selected generated particles"};
  ConfigurableAxis cfgMultiplicityBinning{"cfgMultiplicityBinning", {1000, 0, 10000}, "Binning used for the binning of the number of particles in the event"};
  expressions::Filter zVertexFilter = aod::mccollision::posZ < cfgMaxZVertex && aod::mccollision::posZ > -cfgMaxZVertex;

  void init(InitContext const&)
  {
    mHistManager.add("Multiplicity", "Number of generated particles per MC collision;#bf{#it{N} (Multiplicity)};#bf{#it{N}_{MC collisions}}", HistType::kTH1F, {cfgMultiplicityBinning});
    mHistManager.add("YieldVsMultiplicity", "pT of selected particles in all MC collisions vs number of all generated particles in event;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N} (Multiplicity)};#bf{#it{N}}", HistType::kTH2F, {{200, 0, 20}, cfgMultiplicityBinning});
  }

  void process(soa::Filtered<aod::McCollisions>::iterator const&, aod::McParticles const& mcParticles)
  {
    int nParticles = mcParticles.size();
    mHistManager.fill(HIST("Multiplicity"), nParticles);
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == cfgSelectedParticleCode && std::abs(mcParticle.y()) < cfgRapidityCut)
        mHistManager.fill(HIST("YieldVsMultiplicity"), mcParticle.pt(), nParticles);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<MCGeneratorStudies>(cfgc, TaskName{"mc-generator-studies"})}; }
