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

// Header files.
#include <vector>

// O2 headers.

// O2 Physics headers.
#include "PWGCF/JCorran/Core/FlowJHistManager.h"
#include "CommonConstants/MathConstants.h"

// Namespaces.
using namespace o2;
using namespace o2::framework;

/// \brief Create the histograms in the QA registry.
void FlowJHistManager::createHistQA()
{
  // Security checks for proper use of the method.
  if (!mHistRegistryQA) {
    LOGF(error, "QA histogram registry missing. Quitting...");
    return;
  }

  // Definition of the QA histograms.
  // All the histograms are defined in details for the first centrality
  // class after additional cuts, then cloned for the other classes.
  const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
  mHistRegistryQA->add("Centrality_00-01/After/histCent", "Centrality",
                       HistType::kTH1F, {axisCent}, true);

  const AxisSpec axisMulti{2500, 0., 25000., "N_{tracks}"};
  mHistRegistryQA->add("Centrality_00-01/After/histMulti", "Multiplicity",
                       HistType::kTH1I, {axisMulti}, true);

  const AxisSpec axisZvtx{30, -15., 15., "Z_{vtx} [cm]"};
  mHistRegistryQA->add("Centrality_00-01/After/histZvtx", "Z_{vtx}",
                       HistType::kTH1F, {axisZvtx}, true);

  AxisSpec axisPt = {60, 0., 6., "#it{p}_{T} [GeV/#it{c}]"};
  if (mUseVariablePtBins) {
    std::vector<double> ptBinning = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                     0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                     0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4,
                                     1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6,
                                     2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6.};
    axisPt = {ptBinning, "#it{p}_{T} [GeV/#it{c}]"};
  }
  mHistRegistryQA->add("Centrality_00-01/After/histPt", "#it{p}_{T} (no NUE)",
                       HistType::kTH1F, {axisPt}, true);

  const AxisSpec axisEta = {20, -1., 1., "#eta"};
  mHistRegistryQA->add("Centrality_00-01/After/histEta", "Pseudorapidity",
                       HistType::kTH1F, {axisEta}, true);

  const AxisSpec axisPhi = {100, 0., o2::constants::math::TwoPI, "#varphi"};
  mHistRegistryQA->add("Centrality_00-01/After/histPhi", "Azimuthal angles (no NUA)",
                       HistType::kTH1F, {axisPhi}, true);

  const AxisSpec axisCharge = {2, -2., 2., "Charge"};
  mHistRegistryQA->add("Centrality_00-01/After/histCharge", "Electric charge",
                       HistType::kTH1I, {axisCharge}, true);

  // Additional QA for the full QA task.
  if (mSaveAllQA) {
    // TPC information.
    const AxisSpec axisTPCNcls = {163, -0.5, 162.5, "N_{cls}"};
    mHistRegistryQA->add("Centrality_00-01/After/histTPCNClsFound",
                         "Number of found TPC clusters",
                         HistType::kTH1I, {axisTPCNcls}, true);

    mHistRegistryQA->add("Centrality_00-01/After/histTPCNClsCrossedRows",
                         "Number of crossed TPC rows",
                         HistType::kTH1I, {axisTPCNcls}, true);

    const AxisSpec axisTPCRatio = {20, -0.5, 19.5, "Ratio"};
    mHistRegistryQA->add("Centrality_00-01/After/histTPCCrossedRowsOverFindableCls",
                         "Ratio crossed rows over findable clusters in TPC",
                         HistType::kTH1F, {axisTPCRatio}, true);

    mHistRegistryQA->add("Centrality_00-01/After/histTPCFoundOverFindableCls",
                         "Ratio of found over findable clusters in TPC",
                         HistType::kTH1F, {axisTPCRatio}, true);

    const AxisSpec axisTPCFraction = {30, -0.5, 2.5, "Fraction"};
    mHistRegistryQA->add("Centrality_00-01/After/histTPCFractionSharedCls",
                         "Fraction of shared TPC clusters",
                         HistType::kTH1F, {axisTPCFraction}, true);

    const AxisSpec axisTPCChi = {200, -0.5, 19.5, "#chi^{2} per cl"};
    mHistRegistryQA->add("Centrality_00-01/After/histTPCChi2NCl",
                         "Chi2 per cluster for the TPC track segment",
                         HistType::kTH1F, {axisTPCChi}, true);

    // ITS information.
    const AxisSpec axisITSNcls = {10, -0.5, 9.5, "N_{cls}"};
    mHistRegistryQA->add("Centrality_00-01/After/histITSNCls", "Number of ITS clusters",
                         HistType::kTH1I, {axisITSNcls}, true);

    mHistRegistryQA->add("Centrality_00-01/After/histITSNClsInnerBarrel",
                         "Number of ITS clusters in the Inner Barrel",
                         HistType::kTH1I, {axisITSNcls}, true);

    const AxisSpec axisITSChi = {500, -0.5, 50.5, "#chi^{2} per cl"};
    mHistRegistryQA->add("Centrality_00-01/After/histITSChi2NCl",
                         "Chi2 per cluster for the ITS track segment",
                         HistType::kTH1F, {axisITSChi}, true);

    // DCA information.
    const AxisSpec axisDCAxy = {100, -2.5, 2.5, "DCA_{xy} (cm)"};
    mHistRegistryQA->add("Centrality_00-01/After/histDCAxy", "DCA_{xy} vs #it{p}_{T}",
                         HistType::kTH2F, {axisPt, axisDCAxy}, true);

    const AxisSpec axisDCAz = {110, -5.5, 5.5, "DCA_{z} (cm)"};
    mHistRegistryQA->add("Centrality_00-01/After/histDCAz", "DCA_{z}",
                         HistType::kTH1F, {axisDCAz}, true);

    if (mSaveQABefore) {
      // Clone all the QA for the Before/ distributions.
      mHistRegistryQA->addClone("Centrality_00-01/After/", "Centrality_00-01/Before/");
    }
  }

  if (mObtainNUA) { // TODO: Replace with THnSparse if better with more events.
    mHistRegistryQA->add("Centrality_00-01/After/histZvtxEtaPhi", "Zvtx-eta-phi",
                         HistType::kTH3F, {axisZvtx, axisEta, axisPhi}, true);
  }

  // Add NUE/NUA related histograms for pT and phi only in After/.
  mHistRegistryQA->add("Centrality_00-01/After/histPtCorrected", "#it{p}_{T} (with NUE)",
                       HistType::kTH1F, {axisPt}, true);
  mHistRegistryQA->add("Centrality_00-01/After/histNUEWeights", "NUE weights",
                       HistType::kTH1F, {axisPt}, true);
  mHistRegistryQA->add("Centrality_00-01/After/histPhiCorrected", "Azimuthal angles (with NUA)",
                       HistType::kTH1F, {axisPhi}, true);
  mHistRegistryQA->add("Centrality_00-01/After/histNUAWeights", "NUA weights (projection)",
                       HistType::kTH1F, {axisPhi}, true);

  // Clone the first centrality class into the other classes.
  for (int iBin = 1; iBin < mNcentBins; iBin++) {
    mHistRegistryQA->addClone("Centrality_00-01/", MCentClasses[iBin].data());
  }

  LOGF(info, "QA histograms created.");
}

/// \brief Get the centrality bin value corresponding to the percentile.
/// \param Centrality percentile of the collision.
/// \return Bin for the histograms,...
int FlowJHistManager::getCentBin(float cValue)
{
  const float centClasses[] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 100.};

  for (int i = 0; i < mNcentBins + 1; i++) {
    if (cValue >= centClasses[i]) {
      continue;
    } else {
      return i - 1;
    }
  }

  // We went through all centrality edges without returning at all.
  // --> The measured percentile is larger than the final class we consider.
  return -1;
}
