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

#include "ALICE3/DataModel/OTFMCParticle.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <Framework/AnalysisDataModel.h>

#include <TPDGCode.h>

#include <map>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct Alice3DecayerQA {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    ConfigurableAxis axisCollisionId{"axisCollisionId", {1000, 0, 999}, "CollisionId axis for QA histograms"};
    ConfigurableAxis axisPdgCode{"axisPdgCode", {1000, 0, 999}, "PdgCode axis for QA histograms"};
    ConfigurableAxis axisStatusCode{"axisStatusCode", {1000, 0, 999}, "StatusCode axis for QA histograms"};
    ConfigurableAxis axisFlags{"axisFlags", {10, 0, 9}, "Flags axis for QA histograms"};
    ConfigurableAxis axisMothersIds{"axisMothersIds", {1000, 0, 999}, "MothersIds axis for QA histograms"};
    ConfigurableAxis axisDaughtersIds{"axisDaughtersIds", {1000, 0, 999}, "DaughtersIds axis for QA histograms"};
    ConfigurableAxis axisWeight{"axisWeight", {2, 0, 1}, "Weight axis for QA histograms"};
    ConfigurableAxis axisPos{"axisPos", {1000, 0, 999}, "Position axis for QA histograms"};
    ConfigurableAxis axisPhi{"axisPhi", {720, -360, 360}, "Phi axis for QA histograms"};
    ConfigurableAxis axisEta{"axisEta", {80, -4, 4}, "Eta axis for QA histograms"};
    ConfigurableAxis axisRapidity{"axisRapidity", {80, -4, 4}, "Rapidity axis for QA histograms"};
    ConfigurableAxis axisIsAlive{"axisIsAlive", {2, 0, 1}, "IsAlive axis for QA histograms"};
    ConfigurableAxis axisIsPrimary{"axisIsPrimary", {2, 0, 1}, "IsPrimary axis for QA histograms"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  } axes;

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
    histos.add("DefaultMC/hElPt", "hElPt", kTH1D, {axes.axisPt});
    histos.add("DefaultMC/hMuPt", "hMuPt", kTH1D, {axes.axisPt});
    histos.add("DefaultMC/hPiPt", "hPiPt", kTH1D, {axes.axisPt});
    histos.add("DefaultMC/hKaPt", "hKaPt", kTH1D, {axes.axisPt});
    histos.add("DefaultMC/hPrPt", "hPrPt", kTH1D, {axes.axisPt});

    histos.add("MCWithDau/hElPt", "hElPt", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hMuPt", "hMuPt", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hPiPt", "hPiPt", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hKaPt", "hKaPt", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hPrPt", "hPrPt", kTH1D, {axes.axisPt});

    histos.add("MCWithDau/hCollisionId", "hCollisionId", kTH1D, {axes.axisCollisionId});
    histos.add("MCWithDau/hPdgCode", "hPdgCode", kTH1D, {axes.axisPdgCode});
    histos.add("MCWithDau/hStatusCode", "hStatusCode", kTH1D, {axes.axisStatusCode});
    histos.add("MCWithDau/hFlags", "hFlags", kTH1D, {axes.axisFlags});
    histos.add("MCWithDau/hMothersIds", "hMothersIds", kTH1D, {axes.axisMothersIds});
    histos.add("MCWithDau/hDaughtersIds", "hDaughtersIds", kTH1D, {axes.axisDaughtersIds});
    histos.add("MCWithDau/hWeight", "hWeight", kTH1D, {axes.axisWeight});
    histos.add("MCWithDau/hVx", "hVx", kTH1D, {axes.axisPos});
    histos.add("MCWithDau/hVy", "hVy", kTH1D, {axes.axisPos});
    histos.add("MCWithDau/hVz", "hVz", kTH1D, {axes.axisPos});
    histos.add("MCWithDau/hVt", "hVt", kTH1D, {axes.axisPos});
    histos.add("MCWithDau/hPhi", "hPhi", kTH1D, {axes.axisPhi});
    histos.add("MCWithDau/hEta", "hEta", kTH1D, {axes.axisEta});
    histos.add("MCWithDau/hRapidity", "hRapidity", kTH1D, {axes.axisRapidity});
    histos.add("MCWithDau/hIsAlive", "hIsAlive", kTH1D, {axes.axisIsAlive});
    histos.add("MCWithDau/hIsPrimary", "hIsPrimary", kTH1D, {axes.axisIsPrimary});
    histos.add("MCWithDau/hPx", "hPx", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hPy", "hPy", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hPz", "hPz", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hPt", "hPt", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hP", "hP", kTH1D, {axes.axisPt});
    histos.add("MCWithDau/hE", "hE", kTH1D, {axes.axisPt});
  }

  void processMC(const aod::McParticles&)
  {
    for (const auto& particle : trueEl) {
      histos.fill(HIST("DefaultMC/hElPt"), particle.pt());
    }
    for (const auto& particle : trueMu) {
      histos.fill(HIST("DefaultMC/hMuPt"), particle.pt());
    }
    for (const auto& particle : truePi) {
      histos.fill(HIST("DefaultMC/hPiPt"), particle.pt());
    }
    for (const auto& particle : trueKa) {
      histos.fill(HIST("DefaultMC/hKaPt"), particle.pt());
    }
    for (const auto& particle : truePr) {
      histos.fill(HIST("DefaultMC/hPrPt"), particle.pt());
    }
  }

  void processMCWithDau(const aod::McCollision&, const aod::McParticlesWithDau& particles)
  {
    for (const auto& particle : trueElWithDau) {
      histos.fill(HIST("MCWithDau/hElPt"), particle.pt());
    }
    for (const auto& particle : trueMuWithDau) {
      histos.fill(HIST("MCWithDau/hMuPt"), particle.pt());
    }
    for (const auto& particle : truePiWithDau) {
      histos.fill(HIST("MCWithDau/hPiPt"), particle.pt());
    }
    for (const auto& particle : trueKaWithDau) {
      histos.fill(HIST("MCWithDau/hKaPt"), particle.pt());
    }
    for (const auto& particle : truePrWithDau) {
      histos.fill(HIST("MCWithDau/hPrPt"), particle.pt());
    }

    for (const auto& particle : particles) {
      histos.fill(HIST("MCWithDau/hCollisionId"), particle.mcCollisionId());
      histos.fill(HIST("MCWithDau/hPdgCode"), particle.pdgCode());
      histos.fill(HIST("MCWithDau/hStatusCode"), particle.statusCode());
      histos.fill(HIST("MCWithDau/hFlags"), particle.flags());
      histos.fill(HIST("MCWithDau/hWeight"), particle.weight());
      histos.fill(HIST("MCWithDau/hVx"), particle.vx());
      histos.fill(HIST("MCWithDau/hVy"), particle.vy());
      histos.fill(HIST("MCWithDau/hVz"), particle.vz());
      histos.fill(HIST("MCWithDau/hVt"), particle.vt());
      histos.fill(HIST("MCWithDau/hPhi"), particle.phi());
      histos.fill(HIST("MCWithDau/hEta"), particle.eta());
      histos.fill(HIST("MCWithDau/hRapidity"), particle.y());
      histos.fill(HIST("MCWithDau/hIsAlive"), particle.isAlive());
      histos.fill(HIST("MCWithDau/hIsPrimary"), particle.isPrimary());
      histos.fill(HIST("MCWithDau/hPx"), particle.px());
      histos.fill(HIST("MCWithDau/hPy"), particle.py());
      histos.fill(HIST("MCWithDau/hPz"), particle.pz());
      histos.fill(HIST("MCWithDau/hPt"), particle.pt());
      histos.fill(HIST("MCWithDau/hP"), particle.p());
      histos.fill(HIST("MCWithDau/hE"), particle.e());
      for (const auto& motherParticleId : particle.mothersIds()) {
        histos.fill(HIST("MCWithDau/hMothersIds"), motherParticleId);
      }
      for (const auto& dauParticleId : particle.mothersIds()) {
        histos.fill(HIST("MCWithDau/hDaughtersIds"), dauParticleId);
      }
    }
  }

  PROCESS_SWITCH(Alice3DecayerQA, processMC, "fill MC-only histograms", false);
  PROCESS_SWITCH(Alice3DecayerQA, processMCWithDau, "fill MC-with-dau histograms", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3DecayerQA>(ctx)};
}
