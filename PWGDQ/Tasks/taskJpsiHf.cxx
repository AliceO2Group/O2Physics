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

/// \file taskJPsiHf.cxx
/// \brief Task for the analysis of J/psi - open HF associate production
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// pT differential BDT cuts
namespace bdtcuts
{
static constexpr int nBinsPt = 14;
constexpr float binsPt[nBinsPt + 1] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 16.0f, 24.0f, 36.0f, 1000.0f};
constexpr float bdtCuts[nBinsPt][3] = {{1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}, {1., 0., 0.}};
static const std::vector<std::string> labelsPt{};
static const std::vector<std::string> labelsCutsBdt = {"BDT background", "BDT prompt", "BDT nonprompt"};
} // namespace bdtcuts

// Declarations of various short names
using MyRedEvents = aod::RedJpDmColls;
using MyRedPairCandidatesSelected = aod::RedJpDmDileptons;
using MyRedD0CandidatesSelected = soa::Join<aod::RedJpDmDmesons, aod::RedJpDmD0Masss, aod::RedJpDmDmesBdts, aod::RedJpDmDmDau0s, aod::RedJpDmDmDau1s>;

struct taskJPsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedDleptDmesAll> redDileptDimesAll;

  // HF configurables
  Configurable<float> massHfCandMin{"massHfCandMin", 1.5f, "minimum HF mass"};
  Configurable<float> massHfCandMax{"massHfCandMax", 2.3f, "maximum HF mass"};
  Configurable<std::vector<float>> binsPtDmesForBdt{"binsPtDmesForBdt", std::vector<float>{bdtcuts::binsPt, bdtcuts::binsPt + bdtcuts::nBinsPt + 1}, "pT bin limits for BDT cuts"};
  Configurable<LabeledArray<float>> cutsDmesBdt{"cutsDmesBdt", {bdtcuts::bdtCuts[0], bdtcuts::nBinsPt, 3, bdtcuts::labelsPt, bdtcuts::labelsCutsBdt}, "D-meson BDT selections per pT bin"};

  // DQ configurables
  Configurable<float> massDileptonCandMin{"massDileptonCandMin", 1.f, "minimum dilepton mass"};
  Configurable<float> massDileptonCandMax{"massDileptonCandMax", 5.f, "maximum dilepton mass"};

  // Preslices for unsorted indexes
  PresliceUnsorted<MyRedPairCandidatesSelected> perCollisionDilepton = aod::jpsidmescorr::redJpDmCollId;
  PresliceUnsorted<MyRedD0CandidatesSelected> perCollisionDmeson = aod::jpsidmescorr::redJpDmCollId;

  // histogram for normalisation
  std::shared_ptr<TH1> hCollisions;
  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    hCollisions = registry.add<TH1>("hCollisions", ";;entries", HistType::kTH1F, {{2, -0.5, 1.5}});
    hCollisions->GetXaxis()->SetBinLabel(1, "all collisions");
    hCollisions->GetXaxis()->SetBinLabel(2, "collisions with pairs");
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <typename TEvent, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TEvent const& /*event*/, TDqTrack const& dileptons, THfTrack const& dmesons)
  {
    float ptDilepton = -999;
    float ptDmeson = -999;
    float rapDilepton = -999;
    float rapDmeson = -999;
    float phiDilepton = -999;
    float phiDmeson = -999;
    float deltaRap = -999;
    float deltaPhi = -999;

    for (auto const& dilepton : dileptons) {
      ptDilepton = RecoDecay::pt(dilepton.px(), dilepton.py());
      rapDilepton = RecoDecay::y(std::array{dilepton.px(), dilepton.py(), dilepton.pz()}, constants::physics::MassJPsi);
      phiDilepton = RecoDecay::phi(dilepton.px(), dilepton.py());

      for (auto const& dmeson : dmesons) {
        ptDmeson = RecoDecay::pt(dmeson.px(), dmeson.py());
        phiDmeson = RecoDecay::phi(dmeson.px(), dmeson.py());
        float absDeltaPhiRaw = std::abs(phiDilepton - phiDmeson);
        deltaPhi = (absDeltaPhiRaw < o2::constants::math::PI) ? absDeltaPhiRaw : o2::constants::math::TwoPI - absDeltaPhiRaw;

        auto ptBinDmesForBdt = findBin(binsPtDmesForBdt, ptDmeson);
        if (ptBinDmesForBdt == -1) {
          continue;
        }

        auto minItsClsDmesDau = (dmeson.numItsClsDmesProng0() < dmeson.numItsClsDmesProng1()) ? dmeson.numItsClsDmesProng0() : dmeson.numItsClsDmesProng1();
        auto minTpcCrossRowsDmesDau = (dmeson.numTpcCrossedRowsDmesProng0() < dmeson.numTpcCrossedRowsDmesProng1()) ? dmeson.numTpcCrossedRowsDmesProng0() : dmeson.numTpcCrossedRowsDmesProng1();
        auto minPtDmesDau = (dmeson.ptDmesProng0() < dmeson.ptDmesProng1()) ? dmeson.ptDmesProng0() : dmeson.ptDmesProng1();
        auto minAbsEtaDmesDau = (std::abs(dmeson.etaDmesProng0()) < std::abs(dmeson.etaDmesProng1())) ? std::abs(dmeson.etaDmesProng0()) : std::abs(dmeson.etaDmesProng1());

        if (dmeson.massD0() > 0) {
          rapDmeson = RecoDecay::y(std::array{dmeson.px(), dmeson.py(), dmeson.pz()}, constants::physics::MassD0);
          deltaRap = rapDilepton - rapDmeson;
          auto bdtBkg = dmeson.bdtBkgMassHypo0();
          auto bdtPrompt = dmeson.bdtPromptMassHypo0();
          auto bdtNonPrompt = dmeson.bdtNonpromptMassHypo0();
          if ((dilepton.mass() > massDileptonCandMin && dilepton.mass() < massDileptonCandMax) && (dmeson.massD0() > massHfCandMin && dmeson.massD0() < massHfCandMax && bdtBkg < cutsDmesBdt->get(ptBinDmesForBdt, "BDT background") && bdtPrompt > cutsDmesBdt->get(ptBinDmesForBdt, "BDT prompt") && bdtNonPrompt > cutsDmesBdt->get(ptBinDmesForBdt, "BDT nonprompt"))) {
            redDileptDimesAll(dilepton.mass(), dmeson.massD0(), ptDilepton, ptDmeson, rapDilepton, rapDmeson, phiDilepton, phiDmeson, deltaRap, deltaPhi, bdtBkg, bdtPrompt, bdtNonPrompt, minItsClsDmesDau, minTpcCrossRowsDmesDau, minPtDmesDau, minAbsEtaDmesDau);
          }
        }
        if (dmeson.massD0bar() > 0) {
          rapDmeson = RecoDecay::y(std::array{dmeson.px(), dmeson.py(), dmeson.pz()}, constants::physics::MassD0);
          deltaRap = rapDilepton - rapDmeson;
          auto bdtBkg = dmeson.bdtBkgMassHypo1();
          auto bdtPrompt = dmeson.bdtPromptMassHypo1();
          auto bdtNonPrompt = dmeson.bdtNonpromptMassHypo1();
          if ((dilepton.mass() > massDileptonCandMin && dilepton.mass() < massDileptonCandMax) && (dmeson.massD0bar() > massHfCandMin && dmeson.massD0bar() < massHfCandMax && bdtBkg < cutsDmesBdt->get(ptBinDmesForBdt, "BDT background") && bdtPrompt > cutsDmesBdt->get(ptBinDmesForBdt, "BDT prompt") && bdtNonPrompt > cutsDmesBdt->get(ptBinDmesForBdt, "BDT nonprompt"))) {
            redDileptDimesAll(dilepton.mass(), dmeson.massD0bar(), ptDilepton, ptDmeson, rapDilepton, rapDmeson, phiDilepton, phiDmeson, deltaRap, deltaPhi, bdtBkg, bdtPrompt, bdtNonPrompt, minItsClsDmesDau, minTpcCrossRowsDmesDau, minPtDmesDau, minAbsEtaDmesDau);
          }
        }
      }
    }
  }

  void processRedJspiD0(MyRedEvents const& events, MyRedPairCandidatesSelected const& dileptons, MyRedD0CandidatesSelected const& dmesons)
  {
    // Fill the column of collisions with pairs
    for (auto& event : events) {
      hCollisions->Fill(1.f);
      auto groupedDileptonCandidates = dileptons.sliceBy(perCollisionDilepton, event.index());
      auto groupedDmesonCandidates = dmesons.sliceBy(perCollisionDmeson, event.index());
      runDileptonDmeson(event, groupedDileptonCandidates, groupedDmesonCandidates);
    }
  }

  void processNormCounter(RedJpDmColCounts const& normCounters)
  {
    // Fill the column with all collisions
    for (const auto& normCounter : normCounters) {
      hCollisions->Fill(0.f, static_cast<float>(normCounter.numColls()));
    }
  }

  PROCESS_SWITCH(taskJPsiHf, processRedJspiD0, "Process J/psi - D0", true);
  PROCESS_SWITCH(taskJPsiHf, processNormCounter, "Process normalization counter", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskJPsiHf>(cfgc)};
}
