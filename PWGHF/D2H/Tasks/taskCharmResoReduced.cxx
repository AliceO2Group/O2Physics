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

/// \file taskCharmReso.cxx
/// \brief Charmed Resonances analysis task
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  Ds2StarToDplusK0s,
  XcToDplusLambda,
  LambdaDminus
};

struct HfTaskReso{
// Configurables
ConfigurableAxis configAxisPt{"configAxisPt", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT"};
ConfigurableAxis configAxisPtProng0{"configAxisPtProng0", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT of D daughter"};
ConfigurableAxis configAxisPtProng1{"configAxisPtProng1", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT of V0 daughter"};
ConfigurableAxis configAxisInvMassReso{"configAxisInvMassReso", {225, 2.35, 2.80}, "axis for Resonance invariant mass"};
ConfigurableAxis configAxisInvMassProng0{"configAxisInvMassProng0", {175, 1.70, 2.05}, "axis for Prong0 invariant mass"};
ConfigurableAxis configAxisInvMassProng1{"configAxisInvMassProng1", {80, 0.46, 0.54}, "axis for Prong1 invariant mass"};
ConfigurableAxis configAxisCosThetaStar{"configAxisCosThetaStar", {40, -1, 1}, "axis for V0 costheta star"};
ConfigurableAxis configAxisBkgBdtScore{"configAxisBkgBdtScore", {100, 0, 1}, "axis for D meson Bkg BDT score"};
ConfigurableAxis configAxisNonPromptBdtScore{"configAxisNonPromptBdtScore", {100, 0, 1}, "axis for D meson Non Prompt BDT score"};

using reducedResoWithMl = soa::Join<aod::HfCandCharmReso, aod::HfCharmResoMLs>;
Preslice<aod::HfCandCharmReso> candsResoPerCollision = aod::hf_track_index_reduced::hfRedCollisionId;

// Histogram Registry
HistogramRegistry registry;

//init
void init(InitContext&)
{ 
    const AxisSpec axisPt{configAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtProng0{configAxisPtProng0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtProng1{configAxisPtProng1, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisInvMassReso{configAxisInvMassReso, "inv. mass (D{V}_{0}) (GeV/#it{c}^{2})"};
    const AxisSpec axisInvMassProng0{configAxisInvMassProng0, "inv. mass (D) (GeV/#it{c}^{2})"};
    const AxisSpec axisInvMassProng1{configAxisInvMassProng1, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
    const AxisSpec axisCosThetaStar{configAxisCosThetaStar,"cos(#vartheta*)"};
    const AxisSpec axisBkgBdtScore{configAxisBkgBdtScore,"Bkg BDT Score"};
    const AxisSpec axisNonPromptBdtScore{configAxisNonPromptBdtScore,"Non Prompt BDT Score"};
    registry.add("hMass", "Charm resonance candidates inv. mass", {HistType::kTH1F, {axisInvMassReso}});
    registry.add("hMassProng0", "D daughters inv. mass", {HistType::kTH1F, {axisInvMassProng0}});
    registry.add("hMassProng1", "V0 daughter inv. mass", {HistType::kTH1F, {axisInvMassProng1}});
    registry.add("hPt", "Charm resonance candidates pT", {HistType::kTH1F, {axisPt}});
    registry.add("hPtProng0", "D daughters pT", {HistType::kTH1F, {axisPtProng0}});
    registry.add("hPtProng1", "V0 daughter pT", {HistType::kTH1F, {axisPtProng1}});
    registry.add("hSparseMl", "THn for production studies with cosThStar and BDT scores", HistType::kTHnSparseF,{axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar, axisBkgBdtScore, axisNonPromptBdtScore});
}

// Fill histograms
/// \param channel is the decay channel of the Resonance
/// \param candidate is the candidate table
template <DecayChannel channel, typename Cand>
void fillHisto (const Cand& candidate){
    registry.fill(HIST("hMass"), candidate.invMass());
    registry.fill(HIST("hMassProng0"), candidate.invMassProng0());
    registry.fill(HIST("hMassProng1"), candidate.invMassProng1());
    registry.fill(HIST("hPt"), candidate.pt());
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    float cosThetaStar {0.};
    switch (channel) {
          case DecayChannel::Ds1ToDstarK0s:
            cosThetaStar = candidate.cosThetaStarDs1();
            break;
          case DecayChannel::Ds2StarToDplusK0s:
            cosThetaStar = candidate.cosThetaStarDs2Star();
            break;
          default:
            cosThetaStar = candidate.cosThetaStarXiC3055();
            break;
    }
    registry.fill(HIST("hSparseMl"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.invMass(), candidate.invMassProng0(), candidate.invMassProng1(), cosThetaStar, candidate.mlScoreBkgProng0(), candidate.mlScoreNonpromptProng0());
} // fillHisto

// process functions
void processData(reducedResoWithMl const& candidates)
{
for (const auto& cand : candidates) {
  fillHisto<DecayChannel::Ds2StarToDplusK0s>(cand); 
}
}
PROCESS_SWITCH(HfTaskReso, processData, "Process data", true);

}; //struct HfTaskReso
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskReso>(cfgc)};
}