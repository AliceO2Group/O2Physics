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

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <vector>

using namespace o2;
using namespace o2::framework;

struct Alice3DecayerQA {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  SliceCache cache;

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
    ConfigurableAxis axisRadiusLog{"axisRadiusLog", {VARIABLE_WIDTH, 0.0f, 0.01f, 0.0104713f, 0.0109648f, 0.0114815f, 0.0120226f, 0.0125893f, 0.0131826f, 0.0138038f, 0.0144544f, 0.0151356f, 0.0158489f, 0.0165959f, 0.017378f, 0.018197f, 0.0190546f, 0.0199526f, 0.020893f, 0.0218776f, 0.0229087f, 0.0239883f, 0.0251189f, 0.0263027f, 0.0275423f, 0.0288403f, 0.0301995f, 0.0316228f, 0.0331131f, 0.0346737f, 0.0363078f, 0.0380189f, 0.0398107f, 0.0416869f, 0.0436516f, 0.0457088f, 0.047863f, 0.0501187f, 0.0524807f, 0.0549541f, 0.057544f, 0.060256f, 0.0630957f, 0.0660693f, 0.0691831f, 0.0724436f, 0.0758578f, 0.0794328f, 0.0831764f, 0.0870964f, 0.0912011f, 0.0954993f, 0.1f, 0.104713f, 0.109648f, 0.114815f, 0.120226f, 0.125893f, 0.131826f, 0.138038f, 0.144544f, 0.151356f, 0.158489f, 0.165959f, 0.17378f, 0.18197f, 0.190546f, 0.199526f, 0.20893f, 0.218776f, 0.229087f, 0.239883f, 0.251189f, 0.263027f, 0.275423f, 0.288403f, 0.301995f, 0.316228f, 0.331131f, 0.346737f, 0.363078f, 0.380189f, 0.398107f, 0.416869f, 0.436516f, 0.457088f, 0.47863f, 0.501187f, 0.524807f, 0.549541f, 0.57544f, 0.60256f, 0.630957f, 0.660693f, 0.691831f, 0.724436f, 0.758578f, 0.794328f, 0.831764f, 0.870964f, 0.912011f, 0.954993f, 1.0f, 1.04713f, 1.09648f, 1.14815f, 1.20226f, 1.25893f, 1.31826f, 1.38038f, 1.44544f, 1.51356f, 1.58489f, 1.65959f, 1.7378f, 1.8197f, 1.90546f, 1.99526f, 2.0893f, 2.18776f, 2.29087f, 2.39883f, 2.51189f, 2.63027f, 2.75423f, 2.88403f, 3.01995f, 3.16228f, 3.31131f, 3.46737f, 3.63078f, 3.80189f, 3.98107f, 4.16869f, 4.36516f, 4.57088f, 4.7863f, 5.01187f, 5.24807f, 5.49541f, 5.7544f, 6.0256f, 6.30957f, 6.60693f, 6.91831f, 7.24436f, 7.58578f, 7.94328f, 8.31764f, 8.70964f, 9.12011f, 9.54993f, 10.0f, 10.4713f, 10.9648f, 11.4815f, 12.0226f, 12.5893f, 13.1826f, 13.8038f, 14.4544f, 15.1356f, 15.8489f, 16.5959f, 17.378f, 18.197f, 19.0546f, 19.9526f, 20.893f, 21.8776f, 22.9087f, 23.9883f, 25.1189f, 26.3027f, 27.5423f, 28.8403f, 30.1995f, 31.6228f, 33.1131f, 34.6737f, 36.3078f, 38.0189f, 39.8107f, 41.6869f, 43.6516f, 45.7088f, 47.863f, 50.1187f, 52.4807f, 54.9541f, 57.544f, 60.256f, 63.0957f, 66.0693f, 69.1831f, 72.4436f, 75.8578f, 79.4328f, 83.1764f, 87.0964f, 91.2011f, 95.4993f, 100.0f}, "Radial axis"};
    ConfigurableAxis axisPtLog{"axisPtLog", {VARIABLE_WIDTH, 0.0f, 0.001f, 0.00104713f, 0.00109648f, 0.00114815f, 0.00120226f, 0.00125893f, 0.00131826f, 0.00138038f, 0.00144544f, 0.00151356f, 0.00158489f, 0.00165959f, 0.0017378f, 0.0018197f, 0.00190546f, 0.00199526f, 0.0020893f, 0.00218776f, 0.00229087f, 0.00239883f, 0.00251189f, 0.00263027f, 0.00275423f, 0.00288403f, 0.00301995f, 0.00316228f, 0.00331131f, 0.00346737f, 0.00363078f, 0.00380189f, 0.00398107f, 0.00416869f, 0.00436516f, 0.00457088f, 0.0047863f, 0.00501187f, 0.00524807f, 0.00549541f, 0.0057544f, 0.0060256f, 0.00630957f, 0.00660693f, 0.00691831f, 0.00724436f, 0.00758578f, 0.00794328f, 0.00831764f, 0.00870964f, 0.00912011f, 0.00954993f, 0.01f, 0.0104713f, 0.0109648f, 0.0114815f, 0.0120226f, 0.0125893f, 0.0131826f, 0.0138038f, 0.0144544f, 0.0151356f, 0.0158489f, 0.0165959f, 0.017378f, 0.018197f, 0.0190546f, 0.0199526f, 0.020893f, 0.0218776f, 0.0229087f, 0.0239883f, 0.0251189f, 0.0263027f, 0.0275423f, 0.0288403f, 0.0301995f, 0.0316228f, 0.0331131f, 0.0346737f, 0.0363078f, 0.0380189f, 0.0398107f, 0.0416869f, 0.0436516f, 0.0457088f, 0.047863f, 0.0501187f, 0.0524807f, 0.0549541f, 0.057544f, 0.060256f, 0.0630957f, 0.0660693f, 0.0691831f, 0.0724436f, 0.0758578f, 0.0794328f, 0.0831764f, 0.0870964f, 0.0912011f, 0.0954993f, 0.1f, 0.104713f, 0.109648f, 0.114815f, 0.120226f, 0.125893f, 0.131826f, 0.138038f, 0.144544f, 0.151356f, 0.158489f, 0.165959f, 0.17378f, 0.18197f, 0.190546f, 0.199526f, 0.20893f, 0.218776f, 0.229087f, 0.239883f, 0.251189f, 0.263027f, 0.275423f, 0.288403f, 0.301995f, 0.316228f, 0.331131f, 0.346737f, 0.363078f, 0.380189f, 0.398107f, 0.416869f, 0.436516f, 0.457088f, 0.47863f, 0.501187f, 0.524807f, 0.549541f, 0.57544f, 0.60256f, 0.630957f, 0.660693f, 0.691831f, 0.724436f, 0.758578f, 0.794328f, 0.831764f, 0.870964f, 0.912011f, 0.954993f, 1.0f, 1.04713f, 1.09648f, 1.14815f, 1.20226f, 1.25893f, 1.31826f, 1.38038f, 1.44544f, 1.51356f, 1.58489f, 1.65959f, 1.7378f, 1.8197f, 1.90546f, 1.99526f, 2.0893f, 2.18776f, 2.29087f, 2.39883f, 2.51189f, 2.63027f, 2.75423f, 2.88403f, 3.01995f, 3.16228f, 3.31131f, 3.46737f, 3.63078f, 3.80189f, 3.98107f, 4.16869f, 4.36516f, 4.57088f, 4.7863f, 5.01187f, 5.24807f, 5.49541f, 5.7544f, 6.0256f, 6.30957f, 6.60693f, 6.91831f, 7.24436f, 7.58578f, 7.94328f, 8.31764f, 8.70964f, 9.12011f, 9.54993f, 10.0f}, "pt axis for QA histograms"};
  } axes;

  Partition<aod::McPartWithDaus> trueElectrons = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kElectron);
  Partition<aod::McPartWithDaus> trueMuons = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kMuonMinus);
  Partition<aod::McPartWithDaus> truePions = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kPiPlus);
  Partition<aod::McPartWithDaus> trueKaons = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kKMinus);
  Partition<aod::McPartWithDaus> trueProtons = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kProton);
  Partition<aod::McPartWithDaus> trueK0Short = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kK0Short);
  Partition<aod::McPartWithDaus> trueLambdas = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kLambda0);
  Partition<aod::McPartWithDaus> trueXiMinus = aod::mcparticle::pdgCode == static_cast<int>(PDG_t::kXiMinus);

  template <typename TParticle>
  float radius(const TParticle& particle) const
  {
    return std::hypot(particle.vx(), particle.vy());
  }

  void init(o2::framework::InitContext&)
  {
    // QA with Table entries
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

    // QA with daughters from Decayer
    histos.add("K0S/hGeneratedPt", "hGeneratedPt", kTH1D, {axes.axisPt});
    histos.add("K0S/hHasDecayed", "hHasDecayed", kTH1D, {{2, -0.5, 1.5}});
    histos.add("K0S/hPosDauDecayRadius", "hPosDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("K0S/hNegDauDecayRadius", "hNegDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("Lambda/hGeneratedPt", "hGeneratedPt", kTH1D, {axes.axisPt});
    histos.add("Lambda/hHasDecayed", "hHasDecayed", kTH1D, {{2, -0.5, 1.5}});
    histos.add("Lambda/hPosDauDecayRadius", "hPosDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("Lambda/hNegDauDecayRadius", "hNegDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("XiMinus/hGeneratedPt", "hGeneratedPt", kTH1D, {axes.axisPt});
    histos.add("XiMinus/hHasDecayed", "hHasDecayed", kTH1D, {{2, -0.5, 1.5}});
    histos.add("XiMinus/hBachDauDecayRadius", "hBachDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("XiMinus/hV0DauDecayRadius", "hV0DauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("XiMinus/hPosDauDecayRadius", "hPosDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
    histos.add("XiMinus/hNegDauDecayRadius", "hNegDauDecayRadius", kTH2D, {axes.axisRadiusLog, axes.axisPtLog});
  }

  void process(const aod::McCollision& collision, const aod::McPartWithDaus& particles)
  {
    // Group with collision
    auto trueElectronsGrouped = trueElectrons->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueMuonsGrouped = trueMuons->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto truePionsGrouped = truePions->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueKaonsGrouped = trueKaons->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueProtonsGrouped = trueProtons->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueK0ShortGrouped = trueK0Short->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueLambdasGrouped = trueLambdas->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);
    auto trueXiMinusGrouped = trueXiMinus->sliceByCached(aod::mcparticle::mcCollisionId, collision.globalIndex(), cache);

    for (const auto& particle : trueElectronsGrouped) {
      histos.fill(HIST("MCWithDau/hElPt"), particle.pt());
    }
    for (const auto& particle : trueMuonsGrouped) {
      histos.fill(HIST("MCWithDau/hMuPt"), particle.pt());
    }
    for (const auto& particle : truePionsGrouped) {
      histos.fill(HIST("MCWithDau/hPiPt"), particle.pt());
    }
    for (const auto& particle : trueKaonsGrouped) {
      histos.fill(HIST("MCWithDau/hKaPt"), particle.pt());
    }
    for (const auto& particle : trueProtonsGrouped) {
      histos.fill(HIST("MCWithDau/hPrPt"), particle.pt());
    }
    for (const auto& particle : trueK0ShortGrouped) {
      histos.fill(HIST("K0S/hGeneratedPt"), particle.pt());
      if (particle.has_daughters()) {
        histos.fill(HIST("K0S/hHasDecayed"), 0);
        auto daughters = particle.daughtersIds();
        if (daughters.size() == 2) {
          auto dau0 = particles.rawIteratorAt(daughters.front());
          auto dau1 = particles.rawIteratorAt(daughters.back());

          // K0S -> pi+ pi-
          const bool k0sDecay = (dau0.pdgCode() == PDG_t::kPiPlus && dau1.pdgCode() == PDG_t::kPiMinus) ||
                                (dau0.pdgCode() == PDG_t::kPiMinus && dau1.pdgCode() == PDG_t::kPiPlus);
          if (k0sDecay) {
            auto& positive = dau0.pdgCode() == PDG_t::kPiPlus ? dau0 : dau1;
            auto& negative = dau0.pdgCode() == PDG_t::kPiPlus ? dau1 : dau0;
            histos.fill(HIST("K0S/hPosDauDecayRadius"), radius(positive), positive.pt());
            histos.fill(HIST("K0S/hNegDauDecayRadius"), radius(negative), negative.pt());
          }
        }
      } else {
        histos.fill(HIST("K0S/hHasDecayed"), 1);
      }
    }
    for (const auto& particle : trueLambdasGrouped) {
      histos.fill(HIST("Lambda/hGeneratedPt"), particle.pt());
      if (particle.has_daughters()) {
        histos.fill(HIST("Lambda/hHasDecayed"), 0);
        auto daughters = particle.daughtersIds();
        if (daughters.size() == 2) {
          auto dau0 = particles.rawIteratorAt(daughters[0]);
          auto dau1 = particles.rawIteratorAt(daughters[1]);

          // Lambda -> p pi-
          const bool lambdaDecay = (dau0.pdgCode() == PDG_t::kProton && dau1.pdgCode() == PDG_t::kPiMinus) ||
                                   (dau0.pdgCode() == PDG_t::kPiMinus && dau1.pdgCode() == PDG_t::kProton);
          if (lambdaDecay) {
            auto& positive = dau0.pdgCode() == PDG_t::kProton ? dau0 : dau1;
            auto& negative = dau0.pdgCode() == PDG_t::kProton ? dau1 : dau0;
            histos.fill(HIST("Lambda/hPosDauDecayRadius"), radius(positive), positive.pt());
            histos.fill(HIST("Lambda/hNegDauDecayRadius"), radius(negative), negative.pt());
          }
        }
      } else {
        histos.fill(HIST("Lambda/hHasDecayed"), 1);
      }
    }
    for (const auto& particle : trueXiMinusGrouped) {
      histos.fill(HIST("XiMinus/hGeneratedPt"), particle.pt());
      if (particle.has_daughters()) {
        histos.fill(HIST("XiMinus/hHasDecayed"), 0);
        auto daughters = particle.daughtersIds();
        if (daughters.size() == 2) {
          auto dau0 = particles.rawIteratorAt(daughters.front());
          auto dau1 = particles.rawIteratorAt(daughters.back());

          // Xi- -> Lambda pi-
          const bool xiDecay = (dau0.pdgCode() == PDG_t::kLambda0 && dau1.pdgCode() == PDG_t::kPiMinus) ||
                               (dau0.pdgCode() == PDG_t::kPiMinus && dau1.pdgCode() == PDG_t::kLambda0);
          if (xiDecay) {
            auto& v0 = dau0.pdgCode() == PDG_t::kLambda0 ? dau0 : dau1;
            auto& bachelor = dau0.pdgCode() == PDG_t::kLambda0 ? dau1 : dau0;
            histos.fill(HIST("XiMinus/hBachDauDecayRadius"), radius(bachelor), bachelor.pt());
            histos.fill(HIST("XiMinus/hV0DauDecayRadius"), radius(v0), v0.pt());

            // Lambda -> p pi-
            if (v0.has_daughters()) {
              auto v0daughters = v0.daughtersIds();
              if (v0daughters.size() == 2) {
                auto v0dau0 = particles.rawIteratorAt(v0daughters.front());
                auto v0dau1 = particles.rawIteratorAt(v0daughters.back());
                const bool lambdaDecay = (v0dau0.pdgCode() == PDG_t::kProton && v0dau1.pdgCode() == PDG_t::kPiMinus) ||
                                         (v0dau0.pdgCode() == PDG_t::kPiMinus && v0dau1.pdgCode() == PDG_t::kProton);
                if (lambdaDecay) {
                  auto& positive = v0dau0.pdgCode() == PDG_t::kProton ? v0dau0 : v0dau1;
                  auto& negative = v0dau0.pdgCode() == PDG_t::kProton ? v0dau1 : v0dau0;
                  histos.fill(HIST("XiMinus/hPosDauDecayRadius"), radius(positive), positive.pt());
                  histos.fill(HIST("XiMinus/hNegDauDecayRadius"), radius(negative), negative.pt());
                }
              }
            }
          }
        }
      } else {
        histos.fill(HIST("XiMinus/hHasDecayed"), 1);
      }
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
      for (const auto& dauParticleId : particle.daughtersIds()) {
        histos.fill(HIST("MCWithDau/hDaughtersIds"), dauParticleId);
      }
    }
  }

  PROCESS_SWITCH(Alice3DecayerQA, process, "fill MC-with-dau histograms", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& ctx)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3DecayerQA>(ctx)};
}
