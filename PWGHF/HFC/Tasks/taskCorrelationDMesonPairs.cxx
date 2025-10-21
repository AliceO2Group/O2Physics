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
///
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>, IJCLab Orsay

#include "PWGHF/HFC/DataModel/DMesonPairsTables.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <sys/types.h>

#include <Rtypes.h>

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

enum CandidateType {
  SelectedD = 0, // This particle is selected as a D
  SelectedDbar,  // This particle is selected as a Dbar
  TrueD,         // This particle is a true D
  TrueDbar       // This particle is a true Dbar
};

enum D0Type {
  Default = 0, // Default value
  Signal,      // This particle is a signal D meson
  Reflected,   // This particle is a reflected D meson
  Bkg          // This particle is background of D meson
};
} // namespace

struct HfTaskCorrelationDMesonPairs {
  // Configurables to set sparse axes
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {200, 1.6, 2.1}, "Inv-mass bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {25, 0., 24.}, "pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {10, -1, 1}, "Rapidity bins"};
  ConfigurableAxis thnConfigAxisCandType{"thnConfigAxisCandType", {4, -0.5, 3.5}, "Candidate Type"};

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const AxisSpec thnAxisInvMassCand1{thnConfigAxisInvMass, "inv. mass D_{1} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassCand2{thnConfigAxisInvMass, "inv. mass D_{2} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPtCand1{thnConfigAxisPt, "#it{p}_{T}^{D_{1}} (GeV/#it{c})"};
    const AxisSpec thnAxisPtCand2{thnConfigAxisPt, "#it{p}_{T}^{D_{2}} (GeV/#it{c})"};
    const AxisSpec thnAxisYCand1{thnConfigAxisY, "#it{y} D_{1}"};
    const AxisSpec thnAxisYCand2{thnConfigAxisY, "#it{y} D_{2}"};
    const AxisSpec thnAxisCandType1{thnConfigAxisCandType, "Type D_{1}"};
    const AxisSpec thnAxisCandType2{thnConfigAxisCandType, "Type D_{2}"};

    registry.add("hMass2DCorrelationPairsLS", "hMass2DCorrelationPairsLS", HistType::kTHnSparseD, {thnAxisInvMassCand1, thnAxisInvMassCand2, thnAxisPtCand1, thnAxisPtCand2, thnAxisYCand1, thnAxisYCand2, thnAxisCandType1, thnAxisCandType2});
    registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLS"))->Sumw2();
    registry.add("hMass2DCorrelationPairsOS", "hMass2DCorrelationPairsOS", HistType::kTHnSparseD, {thnAxisInvMassCand1, thnAxisInvMassCand2, thnAxisPtCand1, thnAxisPtCand2, thnAxisYCand1, thnAxisYCand2, thnAxisCandType1, thnAxisCandType2});
    registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOS"))->Sumw2();

    registry.add("hMass2DCorrelationPairsLSMcGen", "hMass2DCorrelationPairsLSMcGen", HistType::kTHnSparseD, {thnAxisInvMassCand1, thnAxisInvMassCand2, thnAxisPtCand1, thnAxisPtCand2, thnAxisYCand1, thnAxisYCand2, thnAxisCandType1, thnAxisCandType2});
    registry.get<THnSparse>(HIST("hMass2DCorrelationPairsLSMcGen"))->Sumw2();
    registry.add("hMass2DCorrelationPairsOSMcGen", "hMass2DCorrelationPairsOSMcGen", HistType::kTHnSparseD, {thnAxisInvMassCand1, thnAxisInvMassCand2, thnAxisPtCand1, thnAxisPtCand2, thnAxisYCand1, thnAxisYCand2, thnAxisCandType1, thnAxisCandType2});
    registry.get<THnSparse>(HIST("hMass2DCorrelationPairsOSMcGen"))->Sumw2();
  }

  uint getD0Type(uint const& candidateType)
  {
    if ((TESTBIT(candidateType, SelectedD) && TESTBIT(candidateType, TrueD)) || (TESTBIT(candidateType, SelectedDbar) && TESTBIT(candidateType, TrueDbar))) {
      return Signal;
    }
    if ((TESTBIT(candidateType, SelectedD) && TESTBIT(candidateType, TrueDbar)) || (TESTBIT(candidateType, SelectedDbar) && TESTBIT(candidateType, TrueD))) {
      return Reflected;
    }
    if ((TESTBIT(candidateType, SelectedD) && !(TESTBIT(candidateType, TrueD) && TESTBIT(candidateType, TrueDbar))) ||
        (TESTBIT(candidateType, SelectedDbar) && !(TESTBIT(candidateType, TrueD) && TESTBIT(candidateType, TrueDbar)))) {
      return Bkg;
    }
    return Default;
  }

  void processData(aod::D0Pair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const ptCand1 = pairEntry.ptCand1();
      float const ptCand2 = pairEntry.ptCand2();
      float const massDCand1 = pairEntry.mDCand1();
      float const massDbarCand1 = pairEntry.mDbarCand1();
      float const massDCand2 = pairEntry.mDCand2();
      float const massDbarCand2 = pairEntry.mDbarCand2();
      float const yCand1 = pairEntry.yCand1();
      float const yCand2 = pairEntry.yCand2();
      auto pairType = pairEntry.pairType();
      auto d0Type1 = getD0Type(pairEntry.candidateType1());
      auto d0Type2 = getD0Type(pairEntry.candidateType2());

      if (TESTBIT(pairType, DD)) {
        registry.fill(HIST("hMass2DCorrelationPairsLS"), massDCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2, d0Type1, d0Type2);
      }
      if (TESTBIT(pairType, DbarDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsLS"), massDbarCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2, d0Type1, d0Type2);
      }

      if (TESTBIT(pairType, DDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsOS"), massDCand1, massDbarCand2, ptCand1, ptCand2, yCand1, yCand2, d0Type1, d0Type2);
      }
      if (TESTBIT(pairType, DbarD)) {
        registry.fill(HIST("hMass2DCorrelationPairsOS"), massDbarCand1, massDCand2, ptCand1, ptCand2, yCand1, yCand2, d0Type1, d0Type2);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processData, "Process data mode", true);

  void processMcGen(aod::D0PairMcGen const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double const ptParticle1 = pairEntry.ptCand1();
      double const ptParticle2 = pairEntry.ptCand2();
      float const massDParticle1 = pairEntry.mDCand1();
      float const massDbarParticle1 = pairEntry.mDbarCand1();
      float const massDParticle2 = pairEntry.mDCand2();
      float const massDbarParticle2 = pairEntry.mDbarCand2();
      float const yParticle1 = pairEntry.yCand1();
      float const yParticle2 = pairEntry.yCand2();
      auto pairType = pairEntry.pairType();
      auto d0Type1 = getD0Type(pairEntry.candidateType1());
      auto d0Type2 = getD0Type(pairEntry.candidateType2());

      if (TESTBIT(pairType, DD)) {
        registry.fill(HIST("hMass2DCorrelationPairsLSMcGen"), massDParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2, d0Type1, d0Type2);
      }
      if (TESTBIT(pairType, DbarDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsLSMcGen"), massDbarParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2, d0Type1, d0Type2);
      }
      if (TESTBIT(pairType, DDbar)) {
        registry.fill(HIST("hMass2DCorrelationPairsOSMcGen"), massDParticle1, massDbarParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2, d0Type1, d0Type2);
      }
      if (TESTBIT(pairType, DbarD)) {
        registry.fill(HIST("hMass2DCorrelationPairsOSMcGen"), massDbarParticle1, massDParticle2, ptParticle1, ptParticle2, yParticle1, yParticle2, d0Type1, d0Type2);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDMesonPairs, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDMesonPairs>(cfgc)};
}
