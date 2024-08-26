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
#include "Common/Core/RecoDecay.h"

// #include "PWGHF/Core/HfHelper.h"
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
// Configurables axis for histos
ConfigurableAxis configAxisPt{"configAxisPt", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT"};
ConfigurableAxis configAxisPtProng0{"configAxisPtProng0", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT of D daughter"};
ConfigurableAxis configAxisPtProng1{"configAxisPtProng1", {VARIABLE_WIDTH,0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "axis for pT of V0 daughter"};
ConfigurableAxis configAxisInvMassReso{"configAxisInvMassReso", {225, 2.35, 2.80}, "axis for Resonance invariant mass"};
ConfigurableAxis configAxisInvMassProng0{"configAxisInvMassProng0", {175, 1.70, 2.05}, "axis for Prong0 invariant mass"};
ConfigurableAxis configAxisInvMassProng1{"configAxisInvMassProng1", {80, 0.46, 0.54}, "axis for Prong1 invariant mass"};
ConfigurableAxis configAxisCosThetaStar{"configAxisCosThetaStar", {40, -1, 1}, "axis for V0 costheta star"};
ConfigurableAxis configAxisBkgBdtScore{"configAxisBkgBdtScore", {100, 0, 1}, "axis for D meson Bkg BDT score"};
ConfigurableAxis configAxisNonPromptBdtScore{"configAxisNonPromptBdtScore", {100, 0, 1}, "axis for D meson Non Prompt BDT score"};
// Configurables for ME
Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
Configurable<int> numberEventsToSkip{"numberEventsToSkip", -1, "Number of events to Skip in ME process"};
ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 50., 100., 200.}, "event multiplicity pools (PV contributors for now)"};
ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "z vertex position pools"};
ConfigurableAxis bzPoolBins{"bzPoolBins", {2, -10, 10}, "Bz of collision"};

using reducedResoWithMl = soa::Join<aod::HfCandCharmReso, aod::HfCharmResoMLs>;
using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib, aod::hf_reduced_collision::Bz>;
SliceCache cache;

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
    registry.add("hNPvCont", "Collision number of PV contributors ; N contrib ; entries", {HistType::kTH1F, {{100, 0, 500}}});
    registry.add("hZvert", "Collision Z Vtx ; z PV [cm] ; entries", {HistType::kTH1F, {{120, -12., 12.}}});
    registry.add("hBz", "Collision Bz ; Bz [T] ; entries", {HistType::kTH1F, {{20, -1., 1.}}});
    if (doprocessDs2StarData){
      registry.add("hSparseMl", "THn for production studies with cosThStar and BDT scores", HistType::kTHnSparseF,{axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar, axisBkgBdtScore, axisNonPromptBdtScore});
    }
    else{ 
      registry.add("hSparse", "THn for production studies with cosThStar", HistType::kTHnSparseF,{axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar});
    }
}

// Fill histograms
/// \tparam channel is the decay channel of the Resonance
/// \tparam hasMl is a flag to fill BDT scores
/// \param candidate is a candidate 
/// \param coll is a reduced collision 
template <DecayChannel channel, bool hasMl, typename Cand, typename Coll>
void fillHisto (const Cand& candidate, const Coll& collision){
    // Collision properties
    registry.fill(HIST("hNPvCont"), collision.numContrib());
    registry.fill(HIST("hZvert"), collision.posZ());
    registry.fill(HIST("hBz"), collision.bz());
    // Candidate properties
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
    if constexpr (hasMl){
          registry.fill(HIST("hSparseMl"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.invMass(), candidate.invMassProng0(), candidate.invMassProng1(), cosThetaStar, candidate.mlScoreBkgProng0(), candidate.mlScoreNonpromptProng0());
    }
    else {
          registry.fill(HIST("hSparse"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.invMass(), candidate.invMassProng0(), candidate.invMassProng1(), cosThetaStar, candidate.mlScoreBkgProng0());
    }
} // fillHisto

// Process data
/// \tparam channel is the decay channel of the Resonance
/// \tparam hasMl is a flag to fill BDT scores
/// \param Coll is the reduced collisions table
/// \param Cand is the candidates table
template <DecayChannel channel, bool hasMl, typename Coll, typename Candidates>
void processData ( Coll const&, Candidates const& candidates){
  for (const auto& cand : candidates) {
    auto coll = cand.template hfRedCollision_as<Coll>();
    fillHisto<channel, hasMl>(cand, coll);
  }
}

// Process data with Mixed Event
/// \tparam channel is the decay channel of the Resonance
/// \tparam hasMl is a flag to fill BDT scores
/// \param Coll is the reduced collisions table
/// \param Cand is the candidates table
template <DecayChannel channel, bool hasMl, typename Coll, typename Candidates>
void processDataME ( Coll const& collisions, Candidates const& candidates){
  BinningType corrBinning{{zPoolBins, multPoolBins, bzPoolBins}, true};
  auto candsTuple = std::make_tuple(candidates);
    SameKindPair<Coll, Candidates, BinningType> pairs{corrBinning, numberEventsMixed, numberEventsToSkip, collisions, candsTuple, &cache}; 
    for (auto& [collision1, cands1, collision2, cands2] : pairs) {
      // For each couple of candidate resonances I can make 2 mixed candidates by swithching daughters
      // To each ME I associate the collision of the D daughter
      for (const auto& [cand1, cand2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(cands1, cands2))) {
        // Vectors with daughter momenta
        std::array<Candidates,2> cands = {cand1, cand2};
        std::array<std::array<float, 3>, 2> pVecsDDau = {cand1.pVectorProng0(),cand2.pVectorProng0()};
        std::array<std::array<float, 3>, 2> pVecsV0Dau = {cand1.pVectorProng1(),cand2.pVectorProng1()};
        std::array<float, 2> bdtBkg = {cand1.mlScoreBkgProng0, cand2.mlScoreBkgProng0()};
        std::array<float, 2> bdtNonPrompt = {cand1.mlScoreNonPromptProng0, cand2.mlScoreNonPromptProng0()};
        // Reconstructing mass, pt and cosThetaStar of mixed candidates
        for (int i = 0; i<2; i++){
          int j = (i + 1) % 2; //index on V0 momenta array
          float ptME = RecoDecay::pt(pVecsDDau[i][0] + pVecsV0Dau[j][0], pVecsDDau[i][1] + pVecsV0Dau[j][1]);
          float invMassME = RecoDecay::m(std::array{pVecsDDau[i], pVecsV0Dau[j]}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0});
          float cosThetaStarME = RecoDecay::cosThetaStar(std::array{pVecsDDau[i], pVecsV0Dau[j]}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0}, invMassME, 1);

        }
      }  
    }
  }

// process functions
void processDs2StarData(aod::HfRedCollisions const& collisions, reducedResoWithMl const& candidates)
{
  processData<DecayChannel::Ds2StarToDplusK0s, true>(collisions, candidates);
}
PROCESS_SWITCH(HfTaskReso, processDs2StarData, "Process data", true);

}; //struct HfTaskReso
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskReso>(cfgc)};
}
