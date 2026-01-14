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

/// \file taskCorrelationDplusDplusReduced.cxx
/// \brief Writer of pairs of D mesons candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Valerio DI BELLA <valerio.di.bella@cern.ch>, IPHC Strasbourg
/// Based on the code of Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/HFC/DataModel/ReducedDMesonPairsTables.h"

#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>
#include <Framework/ASoA.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <cstdlib>
#include <Framework/AnalysisTask.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskCorrelationDplusDplusReduced {
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};

  using SelectedCandidates = soa::Filtered<o2::aod::HfCandDpFulls>;
  using SelectedMcParticles = o2::aod::HfCandDpMcPs;

  Filter filterSelectCandidates = aod::full::candidateSelFlag >= selectionFlagDplus;

  HistogramConfigSpec hTH1NCand{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramConfigSpec hTH1NMcRec{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramConfigSpec hTH1NMcGen{HistType::kTH1F, {{7, -0.5, 6.5}}};
  HistogramRegistry registry{
    "registry",
    {{"hNCand", "Number of D candidates per event;N", hTH1NCand},
     {"hNMcRec", "Number of reconstructed Mc D mesons per event;N", hTH1NMcRec},
     {"hNMcGen", "Number of generated Mc D mesons per event;N", hTH1NMcGen}}};
  void init(InitContext const&)
  {
    registry.add("hMassDplus", "D+ candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDplusMatched", "D+ matched candidates;inv. mass (#pi#pi K) (GeV/#it{c}^{2}))", {HistType::kTH1F, {{120, 1.5848, 2.1848}}});
    registry.add("hMassDMesonPair", "D Meson pair candidates;inv. mass (#pi K) (GeV/#it{c}^{2});inv. mass (#pi K) (GeV/#it{c}^{2})", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {120, 1.5848, 2.1848}}});
    registry.add("hDltPhiMcGen", "Azimuthal correlation for D mesons; #Delta#phi", {HistType::kTH1F, {{100, -3.141593, 3.141593}}});
  }

  void processLocalData(o2::aod::HfCandDpFullEvs::iterator const&,
                        SelectedCandidates const& localCandidates)
  {
    registry.fill(HIST("hNCand"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      auto mass1 = cand1.m();
      for (auto cand2 = cand1 + 1; cand2 != localCandidates.end(); ++cand2) {
        auto mass2 = cand2.m();
        registry.fill(HIST("hMassDMesonPair"), mass2, mass1);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalData, "Process local data", true);

  void processLocalDataMcRec(o2::aod::HfCandDpFullEvs::iterator const&,
                             SelectedCandidates const& localCandidates)
  {
    registry.fill(HIST("hNMcRec"), localCandidates.size());

    for (const auto& cand1 : localCandidates) {
      auto mass1 = cand1.m();
      registry.fill(HIST("hMassDplus"), mass1);
      if (std::abs(cand1.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi)
        registry.fill(HIST("hMassDplusMatched"), mass1);
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalDataMcRec, "Process local MC data", false);

  void processLocalDataMcGen(o2::aod::HfCandDpMcEvs::iterator const&,
                             SelectedMcParticles const& localMcParticles)
  {
    registry.fill(HIST("hNMcGen"), localMcParticles.size());

    for (const auto& part1 : localMcParticles) {
      for (auto part2 = part1 + 1; part2 != localMcParticles.end(); ++part2) {
        registry.fill(HIST("hDltPhiMcGen"), part2.phi() - part1.phi());
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusDplusReduced, processLocalDataMcGen, "Process local MC data at the gen level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskCorrelationDplusDplusReduced>(cfgc)};
}
