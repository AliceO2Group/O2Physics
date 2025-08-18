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
/// \file   mcParticlePrediction.cxx
/// \author Sebastian Scheid, s.scheid@cern.ch
/// \brief Task to build the predictions from the models based on the generated particles
///

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct otfParticlePrediction {
  // histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // define configurables
  ConfigurableAxis binsEta{"binsEta", {100, -5, 5}, "Binning of the Eta axis"};
  ConfigurableAxis binsPt{"binsPt", {100, 0, 10}, "Binning of the Pt axis"};

  Configurable<float> maxYParticle{"maxYParticle", 5.f, "Max rapidity of particles considered"};

  // init function
  void init(InitContext&)
  {

    const AxisSpec axisEta{binsEta, "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    histos.add<TH1>("collisions/generated", "collisions", kTH1D, {{2, -0.5, 1.5}});
    histos.add<TH2>("particles/generated/pi0", "pi0", kTH2D, {axisPt, axisEta});
    histos.add<TH2>("particles/generated/eta", "eta", kTH2D, {axisPt, axisEta});
    histos.add<TH2>("particles/generated/etaP", "etaP", kTH2D, {axisPt, axisEta});
    histos.add<TH2>("particles/generated/rho", "rho", kTH2D, {axisPt, axisEta});
    histos.add<TH2>("particles/generated/omega", "omega", kTH2D, {axisPt, axisEta});
    histos.add<TH2>("particles/generated/phi", "phi", kTH2D, {axisPt, axisEta});
  }

  void process(aod::McCollisions const& mcCollisions,
               aod::McParticles const& mcParticles)
  {

    histos.fill(HIST("collisions/generated"), 0, mcCollisions.size());

    for (const auto& particle : mcParticles) {
      auto pdg = std::abs(particle.pdgCode());
      if (std::abs(particle.y()) > maxYParticle) {
        continue;
      }
      // if (!(particle.isPhysicalPrimary())) {
      //   continue;
      // }
      if (pdg < 100) {
        continue;
      }
      if (pdg > 1000) {
        continue;
      }
      switch (pdg) {
        case 111:
          histos.fill(HIST("particles/generated/pi0"), particle.pt(), particle.y());
          break;
        case 221:
          histos.fill(HIST("particles/generated/eta"), particle.pt(), particle.y());
          break;
        case 331:
          histos.fill(HIST("particles/generated/etaP"), particle.pt(), particle.y());
          break;
        case 223:
          histos.fill(HIST("particles/generated/omega"), particle.pt(), particle.y());
          break;
        case 113:
          histos.fill(HIST("particles/generated/rho"), particle.pt(), particle.y());
          break;
        case 333:
          histos.fill(HIST("particles/generated/phi"), particle.pt(), particle.y());
          break;
        default:
          break;
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<otfParticlePrediction>(cfgc)};
}
