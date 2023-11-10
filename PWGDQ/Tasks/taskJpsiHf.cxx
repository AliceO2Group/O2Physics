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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Declarations of various short names
using MyRedEvents = aod::RedJpDmColls;
using MyRedPairCandidatesSelected = aod::RedJpDmDileptons;
using MyRedD0CandidatesSelected = soa::Join<aod::RedJpDmDmesons, aod::RedJpDmD0Masss, aod::RedJpDmDmesBdts>;

struct taskJPsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedDleptDmesAll> redDileptDimesAll;

  // HF configurables
  Configurable<double> massHfCandMin{"massHfCandMin", 1, "minimum HF mass"};
  Configurable<double> massHfCandMax{"massHfCandMax", 5, "maximum HF mass"};
  // DQ configurables
  Configurable<double> massDileptonCandMin{"massDileptonCandMin", 1, "minimum dilepton mass"};
  Configurable<double> massDileptonCandMax{"massDileptonCandMax", 5, "maximum dilepton mass"};

  // Define histograms
  AxisSpec axisPt{100, 0.f, 50.f};
  AxisSpec axisMassDmeson{200, 1.7f, 2.1f}; // TODO: make it dependent on the D-meson species
  AxisSpec axisMassJPsi{300, 2.f, 5.f};
  AxisSpec axisMidY{60, -1.5f, 1.5f};
  AxisSpec axisFwdY{50, -4.5f, -2.0f};
  AxisSpec axisDeltaY{90, 1.f, 5.5f};
  AxisSpec axisPhi{180, 0., 2 * constants::math::PI}; // same for delta phi

  HistogramRegistry registry{"registry",
                             {{"JPsiDmeson/hSparseJPsiDmeson", ";#it{p}_{T}(J/#psi) (GeV/#it{c});#it{p}_{T}(D) (GeV/#it{c};#it{M}(J/#psi) (GeV/#it{c}^{2});#it{M}(D) (GeV/#it{c}^{2});#Delta#it{y};#Delta#varphi", {HistType::kTHnSparseF, {axisPt, axisPt, axisMassJPsi, axisMassDmeson, axisDeltaY, axisPhi}}},
                              {"JPsiDmeson/hMassJPsiWithDmeson", ";#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassJPsi}}},
                              {"JPsiDmeson/hPtJPsiWithDmeson", ";#it{p}_{T}(J/#psi) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}}},
                              {"JPsiDmeson/hRapJPsiWithDmeson", ";#it{y}(J/#psi);counts", {HistType::kTH1F, {axisFwdY}}},
                              {"JPsiDmeson/hPhiJPsiWithDmeson", ";#it{#varphi}(J/#psi);counts", {HistType::kTH1F, {axisPhi}}},
                              {"JPsiDmeson/hMassDmesonWithJPsi", ";#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassDmeson}}},
                              {"JPsiDmeson/hPtDmesonWithJPsi", ";#it{p}_{T}(D) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}}},
                              {"JPsiDmeson/hRapDmesonWithJPsi", ";#it{y}(D);counts", {HistType::kTH1F, {axisMidY}}},
                              {"JPsiDmeson/hPhiDmesonWithJPsi", ";#it{#varphi}(D);counts", {HistType::kTH1F, {axisPhi}}}}};

  void init(o2::framework::InitContext& context)
  {
    registry.add("JPsi/hMassJPsi", ";#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassJPsi}});
    registry.add("JPsi/hPtJPsi", ";#it{p}_{T}(J/#psi) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}});
    registry.add("JPsi/hRapJPsi", ";#it{y}(J/#psi);counts", {HistType::kTH1F, {axisFwdY}});
    registry.add("JPsi/hPhiJPsi", ";#it{#varphi}(J/#psi);counts", {HistType::kTH1F, {axisPhi}});
    registry.add("Dmeson/hMassDmeson", ";#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassDmeson}});
    registry.add("Dmeson/hPtDmeson", ";#it{p}_{T}(D) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}});
    registry.add("Dmeson/hRapDmeson", ";#it{y}(D);counts", {HistType::kTH1F, {axisMidY}});
    registry.add("Dmeson/hPhiDmeson", ";#it{#varphi}(D);counts", {HistType::kTH1F, {axisPhi}});
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <typename TEvent, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TEvent const& event, TDqTrack const& dileptons, THfTrack const& dmesons)
  {
    float massJpsi = -999;
    float massD0 = -999;
    float massD0bar = -999;

    std::cout << "Dileptons = " << dileptons.size() << " ; D-mesons = " << dmesons.size() << std::endl;
    for (auto& dilepton : dileptons) {
      std::cout << "J/psi collision index = " << dilepton.redJpDmCollId() << std::endl;
      massJpsi = dilepton.mass();
    }
    for (auto& dmeson : dmesons) {
      std::cout << "D0 collision index = " << dmeson.redJpDmCollId() << std::endl;
      massD0    = dmeson.massD0();
      massD0bar = dmeson.massD0bar();
      if (massD0 > 0) {
        registry.fill(HIST("Dmeson/hMassDmeson"), massD0);
      }
      if (massD0bar > 0) {
        registry.fill(HIST("Dmeson/hMassDmeson"), massD0bar);
      }
    }
    std::cout << "----------------" << std::endl;
    redDileptDimesAll(massJpsi, massD0, massD0bar);
  }

  void processRedJspiD0(MyRedEvents::iterator const& event, MyRedPairCandidatesSelected const& dileptons, MyRedD0CandidatesSelected const& dmesons)
  {
    runDileptonDmeson(event, dileptons, dmesons);
  }

  PROCESS_SWITCH(taskJPsiHf, processRedJspiD0, "Process J/psi - D0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskJPsiHf>(cfgc)};
}
