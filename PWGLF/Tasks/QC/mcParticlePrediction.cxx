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
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \brief Task to build the predictions from the models based on the generated particles
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/mcParticle.h"
#include "PWGLF/Utils/inelGt.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;

static const std::vector<std::string> parameterNames{"Enable"};
static constexpr int nParameters = 1;
static const int defaultParameters[o2::pwglf::PIDExtended::NIDsTot][nParameters]{{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
bool enabledArray[o2::pwglf::PIDExtended::NIDsTot];

// Histograms
const int FT0A = 0;
const int FT0C = 1;
const int FT0AC = 2;
const int FV0A = 3;
const int FDDA = 4;
const int FDDC = 5;
const int FDDAC = 6;
const int ZNA = 7;
const int ZNC = 8;
// const int ZEM1 = 9;
// const int ZEM2 = 10;
// const int ZPA = 11;
// const int ZPC = 12;
const int nEstimators = 13;

const char* estimatorNames[nEstimators] = {"FT0A",
                                           "FT0C",
                                           "FT0AC",
                                           "FV0A",
                                           "FDDA",
                                           "FDDC",
                                           "FDDAC",
                                           "ZNA",
                                           "ZNC",
                                           "ZEM1",
                                           "ZEM2",
                                           "ZPA",
                                           "ZPC"};

std::array<std::array<std::shared_ptr<TH2>, o2::pwglf::PIDExtended::NIDsTot>, nEstimators> hpt;
std::array<std::array<std::shared_ptr<TH1>, o2::pwglf::PIDExtended::NIDsTot>, nEstimators> hyield;

struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsPt{"binsPt", {100, 0, 10}, "Binning of the Pt axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 1000}, "Binning of the Multiplicity axis"};
  Configurable<LabeledArray<int>> enabledSpecies{"enabledSpecies",
                                                 {defaultParameters[0], o2::pwglf::PIDExtended::NIDsTot, nParameters, o2::pwglf::PIDExtended::arrayNames(), parameterNames},
                                                 "Bethe Bloch parameters"};
  Configurable<bool> selectInelGt0{"selectInelGt0", true, "Select only inelastic events"};
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(o2::framework::InitContext&)
  {
    histos.add("collisions", "collisions", kTH1D, {{10, 0, 10}});
    auto h = histos.add<TH1>("particles", "particles", kTH1D, {{o2::pwglf::PIDExtended::NIDsTot, -0.5, -0.5 + o2::pwglf::PIDExtended::NIDsTot}});
    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, o2::pwglf::PIDExtended::getName(i));
    }
    h = histos.add<TH1>("particlesPrim", "particlesPrim", kTH1D, {{o2::pwglf::PIDExtended::NIDsTot, -0.5, -0.5 + o2::pwglf::PIDExtended::NIDsTot}});
    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, o2::pwglf::PIDExtended::getName(i));
    }

    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisMultiplicity{binsMultiplicity, "Mult"};
    // FT0
    histos.add("multiplicity/FT0A", "FT0A", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/FT0C", "FT0C", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/FT0AC", "FT0AC", kTH1D, {axisMultiplicity});

    // FV0
    histos.add("multiplicity/FV0A", "FV0A", kTH1D, {axisMultiplicity});

    // FDD
    histos.add("multiplicity/FDDA", "FDDA", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/FDDC", "FDDC", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/FDDAC", "FDDAC", kTH1D, {axisMultiplicity});

    // ZDC
    histos.add("multiplicity/ZNA", "ZNA", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/ZNC", "ZNC", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/ZEM1", "ZEM1", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/ZEM2", "ZEM2", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/ZPA", "ZPA", kTH1D, {axisMultiplicity});
    histos.add("multiplicity/ZPC", "ZPC", kTH1D, {axisMultiplicity});

    // ITS
    histos.add("multiplicity/ITS", "ITS", kTH1D, {axisMultiplicity});

    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      if (enabledSpecies->get(o2::pwglf::PIDExtended::getName(i), "Enable") != 1) {
        enabledArray[i] = false;
        continue;
      }
      enabledArray[i] = true;
      for (int j = 0; j < nEstimators; j++) {
        hpt[j][i] = histos.add<TH2>(Form("pt/%s/%s", estimatorNames[j], o2::pwglf::PIDExtended::getName(i)), o2::pwglf::PIDExtended::getName(i), kTH2D, {binsPt, axisMultiplicity});
        hyield[j][i] = histos.add<TH1>(Form("yield/%s/%s", estimatorNames[j], o2::pwglf::PIDExtended::getName(i)), o2::pwglf::PIDExtended::getName(i), kTH1D, {axisMultiplicity});
      }
    }
  }

  int countInAcceptance(const aod::McParticles& mcParticles, const float etamin, const float etamax)
  {
    // static_assert(etamin < etamax, "etamin must be smaller than etamax");
    int counter = 0;
    for (const auto& particle : mcParticles) {
      if (particle.eta() > etamin && particle.eta() < etamax) {
        counter++;
      }
    }
    return counter;
  }

  int countFT0A(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, 3.5f, 4.9f); }
  int countFT0C(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, -3.3f, -2.1f); }
  int countFV0A(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, 2.2f, 5.1f); }
  int countFDDA(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, 4.9f, 6.3f); }
  int countFDDC(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, -7.f, -4.9f); }
  int countZNA(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, 8.8f, 100.f); }
  int countZNC(const aod::McParticles& mcParticles) { return countInAcceptance(mcParticles, -100.f, -8.8f); }

  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisions"), 0.5);
    const int nFT0A = countFT0A(mcParticles);
    const int nFT0C = countFT0C(mcParticles);
    const int nFV0A = countFV0A(mcParticles);
    const int nFDDA = countFDDA(mcParticles);
    const int nFDDC = countFDDC(mcParticles);
    const int nZNA = countZNA(mcParticles);
    const int nZNC = countZNC(mcParticles);

    for (const auto& particle : mcParticles) {
      particle.pdgCode();
      const auto id = o2::pwglf::PIDExtended::pdgToId(particle);
      if (id < 0) {
        continue;
      }
      if (!enabledArray[id]) {
        continue;
      }
      hpt[FT0A][id]->Fill(particle.pt(), nFT0A);
      hpt[FT0C][id]->Fill(particle.pt(), nFT0C);
      hpt[FT0AC][id]->Fill(particle.pt(), nFT0A + nFT0C);

      hpt[FV0A][id]->Fill(particle.pt(), nFV0A);
      hpt[FDDA][id]->Fill(particle.pt(), nFDDA);
      hpt[FDDC][id]->Fill(particle.pt(), nFDDC);
      hpt[FDDAC][id]->Fill(particle.pt(), nFDDA + nFDDC);
      hpt[ZNA][id]->Fill(particle.pt(), nZNA);
      hpt[ZNC][id]->Fill(particle.pt(), nZNC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
