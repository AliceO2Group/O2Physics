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

/// \file taskCorrelationDMesonPairs.cxx
/// \brief D-Dbar analysis task - data-like, MC-reco and MC-kine analyses for D0 and DPlus pairs.
/// \note Temporary code for tests of D0-D0 correlations
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/DMesonPairsTablesTesting.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
enum PairTypeSel {
  DD = 0,   // D0-D0
  DbarDbar, // D0bar-D0bar
  DDbar,
  DbarD
};
} // namespace

// string definitions, used for histogram axis labels
const char stringCorrelationPairs[180] = "D meson pair candidates 2D;inv. mass D_{1} (GeV/#it{c}^{2});inv. mass D_{2} (GeV/#it{c}^{2});#it{p}_{T}^{D_{1}} (GeV/#it{c});#it{p}_{T}^{D_{2}} (GeV/#it{c});#eta D_{1};#eta D_{2};";
const char stringCorrelationPairsFinerBinning[194] = "D meson pair candidates 2D Finer Binning;inv. mass D_{1} (GeV/#it{c}^{2});inv. mass D_{2} (GeV/#it{c}^{2});#it{p}_{T}^{D_{1}} (GeV/#it{c});#it{p}_{T}^{D_{2}} (GeV/#it{c});#eta D_{1};#eta D_{2};";

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 8;
const double ptBinsCorrelations[nPtBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecPtBinsCorrelations = std::vector<double>{ptBinsCorrelations, ptBinsCorrelations + nPtBinsCorrelations + 1};

struct HfTaskCorrelationDMesonPairsTesting {
  // Enable histograms with finer pT and y binning
  Configurable<bool> enableFinerBinning{"enableFinerBinning", false, "Enable histograms with finer pT and y binning"};
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecPtBinsCorrelations}, "pT bin limits for correlation plots"};

  // HistoTypes
  HistogramConfigSpec hTHnMass2DCorrPairs{HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {10, 0., 10.}, {10, 0., 10.}, {10, -1, 1}, {10, -1, 1}}}; // note: axes 3 and 4 (the pT) are updated in the init();
  HistogramConfigSpec hTHnMass2DCorrPairsFinerBinning{HistType::kTHnSparseD, {{200, 1.6, 2.1}, {200, 1.6, 2.1}, {60, 1., 6.}, {60, 1., 6.}, {160, -0.8, 0.8}, {160, -0.8, 0.8}}};

  HistogramRegistry registry{
    "registry",
    {{"hMass2DCorrelationPairsLS", stringCorrelationPairs, hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsOS", stringCorrelationPairs, hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsLSMcGen", stringCorrelationPairs, hTHnMass2DCorrPairs},
     {"hMass2DCorrelationPairsOSMcGen", stringCorrelationPairs, hTHnMass2DCorrPairs}}};

  void init(InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinspTaxis = binsPtCorrelations->size() - 1;
    const double* valuespTaxis = binsPtCorrelations->data();

    if (enableFinerBinning) {
      registry.add("hMass2DCorrelationPairsLSFinerBinning", stringCorrelationPairsFinerBinning, hTHnMass2DCorrPairsFinerBinning);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLSFinerBinning"))->Sumw2();
      registry.add("hMass2DCorrelationPairsOSFinerBinning", stringCorrelationPairsFinerBinning, hTHnMass2DCorrPairsFinerBinning);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOSFinerBinning"))->Sumw2();

      registry.add("hMass2DCorrelationPairsLSMcGenFinerBinning", stringCorrelationPairsFinerBinning, hTHnMass2DCorrPairsFinerBinning);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLSMcGenFinerBinning"))->Sumw2();
      registry.add("hMass2DCorrelationPairsOSMcGenFinerBinning", stringCorrelationPairsFinerBinning, hTHnMass2DCorrPairsFinerBinning);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOSMcGenFinerBinning"))->Sumw2();
    }

    for (int i = 2; i <= 3; i++) {
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLS"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLS"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOS"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOS"))->Sumw2();

      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLSMcGen"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLSMcGen"))->Sumw2();
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOSMcGen"))->GetAxis(i)->Set(nBinspTaxis, valuespTaxis);
      registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOSMcGen"))->Sumw2();
    }
  }

  void processData(aod::D0PairTesting const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float ptCand1 = pairEntry.ptCand1();
      float ptCand2 = pairEntry.ptCand2();
      float massDCand1 = pairEntry.mDCand1();
      float massDbarCand1 = pairEntry.mDbarCand1();
      float massDCand2 = pairEntry.mDCand2();
      float massDbarCand2 = pairEntry.mDbarCand2();
      float yCand1 = pairEntry.yCand1();
      float yCand2 = pairEntry.yCand2();
      auto pairType = pairEntry.pairType();

      if (TESTBIT(pairType, DD)) {
        registry.fill(HIST("hMass2DCorrelationPairsLS"), massDCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsLSFinerBinning"), massDCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      }
      if (TESTBIT(pairType, DbarDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsLS"), massDbarCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsLSFinerBinning"), massDbarCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      }

      if (TESTBIT(pairType, DDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsOS"), massDCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsOSFinerBinning"), massDCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      }
      if (TESTBIT(pairType, DbarD)) {
        registry.fill(HIST("hMass2DCorrelationPairsOS"), massDbarCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsOSFinerBinning"), massDbarCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairsTesting, processData, "Process data mode", true);

  void processMcGen(aod::D0PairMcGenTesting const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double ptParticle1 = pairEntry.ptCand1();
      double ptParticle2 = pairEntry.ptCand2();
      float massDParticle1 = pairEntry.mDCand1();
      float massDbarParticle1 = pairEntry.mDbarCand1();
      float massDParticle2 = pairEntry.mDCand2();
      float massDbarParticle2 = pairEntry.mDbarCand2();
      float yParticle1 = pairEntry.yCand1();
      float yParticle2 = pairEntry.yCand2();
      auto pairType = pairEntry.pairType();

      if (TESTBIT(pairType, DD)) {
        registry.fill(HIST("hMass2DCorrelationPairsLSMcGen"), massDParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsLSMcGenFinerBinning"), massDParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        }
      }
      if (TESTBIT(pairType, DbarDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsLSMcGen"), massDbarParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsLSMcGenFinerBinning"), massDbarParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        }
      }
      if (TESTBIT(pairType, DDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsOSMcGen"), massDParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsOSMcGenFinerBinning"), massDParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        }
      }
      if (TESTBIT(pairType, DbarD)) {
        registry.fill(HIST("hMass2DCorrelationPairsOSMcGen"), massDbarParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        if (enableFinerBinning) {
          registry.fill(HIST("hMass2DCorrelationPairsOSMcGenFinerBinning"), massDbarParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairsTesting, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDMesonPairsTesting>(cfgc)};
}
