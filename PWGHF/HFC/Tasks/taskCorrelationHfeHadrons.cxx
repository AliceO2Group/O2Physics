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

/// \file taskCorrelationHfeHadrons.cxx
/// \brief HFE-Hadrons azimuthal correlations analysis task - data-like, MC-reco and MC-Gen analyses
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "PWGHF/HFC/DataModel/CorrelationTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_correlation_electron_hadron;

struct HfTaskCorrelationHfeHadrons {
  // Configurables
  // Deltaphi binning
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 32, "Bins for #Delta#varphi bins"};

  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {30, -1.8, 1.8}, "#it{#Delta#eta}"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {32, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "#it{#Delta#varphi}"};
  ConfigurableAxis binsPt{"binsPt", {50, 0.0, 50}, "#it{p_{T}}(GeV/#it{c})"};

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta #eta = #eta_{Electron}- #eta_{Hadron}"};
    AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta #varphi = #varphi_{Electron}- #varphi_{Hadron}"};
    AxisSpec axisPt = {binsPt, "#it{p_{T}}(GeV/#it{c})"};

    registry.add("hInclusiveEHCorrel", "Sparse for Delta phi and Delta eta Inclusive Electron with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hLikeSignEHCorrel", "Sparse for Delta phi and Delta eta Like sign Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hUnLikeSignEHCorrel", "Sparse for Delta phi and Delta eta  UnLike sign Electron pair with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMcGenInclusiveEHCorrel", "Sparse for Delta phi and Delta eta for McGen Inclusive Electron  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hMcGenNonHfEHCorrel", "Sparse for Delta phi and Delta eta  for McGen  NonHeavy flavour Electron pair  with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", {HistType::kTHnSparseF, {{axisPt}, {axisPt}, {axisDeltaPhi}, {axisDeltaEta}}});
  }

  // correlation  for electron hadron
  void process(aod::HfEHadronPair const& pairEntries)
  {
    double deltaPhi = -999;
    double deltaEta = -999;
    double ptHadron = -999;
    double ptElectron = -999;

    for (const auto& pairEntry : pairEntries) {

      deltaPhi = pairEntry.deltaPhi();
      deltaEta = pairEntry.deltaEta();
      ptElectron = pairEntry.ptElectron();
      ptHadron = pairEntry.ptHadron();

      registry.fill(HIST("hInclusiveEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
      if (pairEntry.nPairsLS() > 0) {
        for (int i = 0; i < pairEntry.nPairsLS(); ++i) {

          registry.fill(HIST("hLikeSignEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
        }
      }
      if (pairEntry.nPairsUS() > 0) {
        for (int i = 0; i < pairEntry.nPairsLS(); ++i) {

          registry.fill(HIST("hUnlikeSignEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
        }
      }
    }
  }

  PROCESS_SWITCH(HfTaskCorrelationHfeHadrons, process, "Process ", false);

  void processMcGen(aod::HfEHadronMcPair const& mcGenpairEntries)
  {
    double deltaPhi = -999;
    double deltaEta = -999;
    double ptHadron = -999;
    double ptElectron = -999;

    for (const auto& pairEntry : mcGenpairEntries) {

      deltaPhi = pairEntry.deltaPhi();
      deltaEta = pairEntry.deltaEta();
      ptElectron = pairEntry.ptElectron();
      ptHadron = pairEntry.ptHadron();

      registry.fill(HIST("hMcGenInclusiveEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
      if (pairEntry.isNonHfEHCorr() != 0) {

        registry.fill(HIST("hMcGenNonHfEHCorrel"), ptElectron, ptHadron, deltaPhi, deltaEta);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationHfeHadrons, processMcGen, "Process for Mc Gen ", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationHfeHadrons>(cfgc)};
}
