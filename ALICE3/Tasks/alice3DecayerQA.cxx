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
/// \file alice3DecayerQA.cxx
///
/// \brief QA task for otf decayer
///
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch> 
/// \since  Dec 23, 2025
///

#include <map>
#include <vector>

#include <Framework/AnalysisDataModel.h>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ASoAHelpers.h"

#include "ALICE3/DataModel/OTFMCParticle.h"

#include <TPDGCode.h>


using namespace o2;
using namespace o2::framework;

struct Alice3DecayerQA {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.22f, 1.42f}, ""};
  ConfigurableAxis axisRadius{"axisRadius", {1000, 0, 100}, "Radius"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Partition<aod::McParticles> trueEl = aod::mcparticle::pdgCode == static_cast<int>(kElectron);
  Partition<aod::McParticles> trueMu = aod::mcparticle::pdgCode == static_cast<int>(kMuonMinus);
  Partition<aod::McParticles> truePi = aod::mcparticle::pdgCode == static_cast<int>(kPiPlus);
  Partition<aod::McParticles> trueKa = aod::mcparticle::pdgCode == static_cast<int>(kKMinus);
  Partition<aod::McParticles> truePr = aod::mcparticle::pdgCode == static_cast<int>(kProton);
  
  Partition<aod::McParticlesWithDau> trueElWithDau = aod::mcparticle::pdgCode == static_cast<int>(kElectron);
  Partition<aod::McParticlesWithDau> trueMuWithDau = aod::mcparticle::pdgCode == static_cast<int>(kMuonMinus);
  Partition<aod::McParticlesWithDau> truePiWithDau = aod::mcparticle::pdgCode == static_cast<int>(kPiPlus);
  Partition<aod::McParticlesWithDau> trueKaWithDau = aod::mcparticle::pdgCode == static_cast<int>(kKMinus);
  Partition<aod::McParticlesWithDau> truePrWithDau = aod::mcparticle::pdgCode == static_cast<int>(kProton);


  void init(o2::framework::InitContext&)
  {
    histos.add("DefaultMC/hElPt", "hElPt", kTH1D, {axisPt});
    histos.add("DefaultMC/hMuPt", "hMuPt", kTH1D, {axisPt});
    histos.add("DefaultMC/hPiPt", "hPiPt", kTH1D, {axisPt});
    histos.add("DefaultMC/hKaPt", "hKaPt", kTH1D, {axisPt});
    histos.add("DefaultMC/hPrPt", "hPrPt", kTH1D, {axisPt});

    histos.add("MCWithDau/hElPt", "hElPt", kTH1D, {axisPt});
    histos.add("MCWithDau/hMuPt", "hMuPt", kTH1D, {axisPt});
    histos.add("MCWithDau/hPiPt", "hPiPt", kTH1D, {axisPt});
    histos.add("MCWithDau/hKaPt", "hKaPt", kTH1D, {axisPt});
    histos.add("MCWithDau/hPrPt", "hPrPt", kTH1D, {axisPt});

    histos.add("Lambda/hGenLambda", "hGenLambda", kTH2D, {axisRadius, axisPt});
    histos.add("Xi/hGenXi", "hGenXi", kTH2D, {axisRadius, axisPt});
  }

  void processMC(const aod::McParticles&)
  {
    for (auto const& particle : trueEl) {
      histos.fill(HIST("DefaultMC/hElPt"), particle.pt());
    }
    for (auto const& particle : trueMu) {
      histos.fill(HIST("DefaultMC/hMuPt"), particle.pt());
    }
    for (auto const& particle : truePi) {
      histos.fill(HIST("DefaultMC/hPiPt"), particle.pt());
    }
    for (auto const& particle : trueKa) {
      histos.fill(HIST("DefaultMC/hKaPt"), particle.pt());
    }
    for (auto const& particle : truePr) {
      histos.fill(HIST("DefaultMC/hPrPt"), particle.pt());
    }
  }

  void processMCWithDau(const aod::McParticlesWithDau&)
  {
    for (auto const& particle : trueElWithDau) {
      histos.fill(HIST("MCWithDau/hElPt"), particle.pt());
    }
    for (auto const& particle : trueMuWithDau) {
      histos.fill(HIST("MCWithDau/hMuPt"), particle.pt());
    }
    for (auto const& particle : truePiWithDau) {
      histos.fill(HIST("MCWithDau/hPiPt"), particle.pt());
    }
    for (auto const& particle : trueKaWithDau) {
      histos.fill(HIST("MCWithDau/hKaPt"), particle.pt());
    }
    for (auto const& particle : truePrWithDau) {
      histos.fill(HIST("MCWithDau/hPrPt"), particle.pt());
    }
  }

  PROCESS_SWITCH(Alice3DecayerQA, processMC, "fill MC-only histograms", true);
  PROCESS_SWITCH(Alice3DecayerQA, processMCWithDau, "fill MC-only histograms", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3DecayerQA>(ctx)};
}
