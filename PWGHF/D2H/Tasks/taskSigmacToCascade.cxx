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

/// \file taskSigmacToCascade.cxx
/// \brief Task for Σc0,++ → Λc(→ K0sP) + π-,+ analysis
/// \note Here the Lc from the cascade channel is obtained using the task taskLcToK0sP.cxx
/// \author Rutuparna Rath <rrath@cern.ch>, INFN BOLOGNA and GSI Darmstadt
/// In collaboration with Andrea Alici <aalici@cern.ch>, INFN BOLOGNA

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskSigmacToCascade {

  Configurable<bool> enableTHn{"enableTHn", false, "enable the usage of THn for Σc0, Σc++  and Σc0,++"};
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<int> selectionFlagLcbarToK0sP{"selectionFlagLcbarToK0sP", 1, "Selection Flag for Lcbar"};
  Configurable<float> etaCandMax{"etaCandMax", -1., "max. cand. pseudorapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {16, 0, 16}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisDecLengthXY{"thnConfigAxisDecLengthXY", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisCPAXY{"thnConfigAxisCPAXY", {20, 0.8, 1}, ""};
  ConfigurableAxis configAxisMassLambdaC{"configAxisMassLambdaC", {600, 1.98, 2.58}, ""};
  ConfigurableAxis configAxisDeltaMassSigmaC{"configAxisDeltaMassSigmaC", {200, 0.13, 0.23}, ""};
  Configurable<float> nSigmaSoftPi{"pionNSigma", 3., "NSigma TPC selection"};
  ConfigurableAxis configAxisChargeSigmaC{"configAxisChargeSigmaC", {3, 0, 3}, "charge of SigmaC"};

  /// Filter the candidate Λc+ used for the Σc0,++ creation
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP ||
                                    aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  using RecoLc = soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP>>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    // axes
    AxisSpec const axisBinsPt = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisPt = {300, 0.0f, 30.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisEta = {500, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec const axisY = {500, -2.0f, 2.0f, "y"};
    AxisSpec const axisPhi = {100, 0.f, 6.3f, "#it{#phi}"};
    AxisSpec const axisMassCand = {600, 1.98f, 2.58f, "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisd0 = {500, -0.5f, 0.5f, "DCAxy (cm)"};
    AxisSpec const axisd0V0Daughters = {1000, -5.0f, 5.0f, "DCAxy (cm)"};
    AxisSpec const axisV0CPA = {500, 0.98f, 1.0001f, "v0 cos pointing angle"};
    AxisSpec const axisV0Radius = {1000, 0.f, 40.f, "V0 radius (cm)"};
    AxisSpec const axisV0DCADaughters = {200, 0.f, 2.f, "DCA (cm)"};
    AxisSpec const axisMassK0Short = {500, 0.4f, 0.6f, "#it{m}(K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassLambda = {500, 1.0f, 1.2f, "#it{m}(#Lambda) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassGamma = {500, 0.0f, 0.4f, "#it{m}(#gamma) (GeV/#it{c}^{2})"};
    AxisSpec const axisCPACand = {110, -1.1f, 1.1f, "candiate cos pointing angle"};
    AxisSpec const axisDecLength = {200, 0.f, 2.0f, "decay length (cm)"};
    AxisSpec const axisProperLifetime = {100, 0.f, 0.2f, "#it{c#tau} (cm)"};
    AxisSpec const axisProperLifetimeV0 = {1000, 0.f, 80.f, "#it{c#tau} (cm)"};
    AxisSpec const axisNSigma = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec const axisPidP = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisDeltaMassSigmaC{configAxisDeltaMassSigmaC, "#it{M}(pK_{S}^{0}#pi) - #it{M}(pK_{S}^{0}) (GeV/#it{c}^{2})"};

    // data
    ////////////////////////////////////////////////////////////////////////////
    /// Declare histograms related to Sigma_C analysis from LcToK0sP channel///
    ///////////////////////////////////////////////////////////////////////////

    registry.add("Data/hDeltaMassSc0", "#Sigma_{c}^{0} candidates; #it{M}(K0sP #pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); counts;", {HistType::kTH1D, {{200, 0.13, 0.23}}});                                                                  /// Σc0
    registry.add("Data/hDeltaMassScPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(K0sP#pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); counts;", {HistType::kTH1D, {{200, 0.13, 0.23}}});                                                           /// Σc++
    registry.add("Data/hDeltaMassSc0PlusPlus", "#Sigma_{c}^{0, ++} candidates; #it{M}(K0sP#pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); counts;", {HistType::kTH1D, {{200, 0.13, 0.23}}});                                                       /// Σc0,++
    registry.add("Data/hDeltaMassSc0VsPt", "#Sigma_{c}^{0} candidates; #it{M}(K0sP #pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});              /// Σc0
    registry.add("Data/hDeltaMassScPlusPlusVsPt", "#Sigma_{c}^{++} candidates; #it{M}(K0sP #pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});     /// Σc++
    registry.add("Data/hDeltaMassSc0PlusPlusVsPt", "#Sigma_{c}^{0, ++} candidates; #it{M}(K0sP #pi) - #it{M}(K0sP) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}); /// Σc0, ++
    registry.add("Data/hEtaSc0", "#Sigma_{c}^{0} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                                                                                             /// Σc0
    registry.add("Data/hEtaScPlusPlus", "#Sigma_{c}^{++} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                                                                                     /// Σc++
    registry.add("Data/hEtaSc0PlusPlus", "#Sigma_{c}^{0, ++} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                                                                                 /// Σc0,++
    registry.add("Data/hEtaSc0PlusPlusVsPt", "#Sigma_{c}^{0, ++} candidates; #eta; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisEta, axisPt}});                                                                     /// Σc0,++
    registry.add("Data/hYSc0PlusPlusVsPt", "#Sigma_{c}^{0, ++} candidates; y (rapidity); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisY, axisPt}});                                                                 /// Σc0,++
    registry.add("Data/hPhiSc0", "#Sigma_{c}^{0} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                                                                                             /// Σc0
    registry.add("Data/hPhiScPlusPlus", "#Sigma_{c}^{++} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                                                                                     /// Σc++
    registry.add("Data/hPhiSc0PlusPlus", "#Sigma_{c}^{0, ++} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                                                                                 /// Σc0,++

    /// softPions
    registry.add("Data/hPtSoftPiSc0", "#pi^{-} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{-} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); counts;", {HistType::kTH1D, {axisPt}});                 /// Σc0
    registry.add("Data/hPtSoftPiScPlusPlus", "#pi^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); counts;", {HistType::kTH1D, {axisPt}});         /// Σc++
    registry.add("Data/hPtSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0, ++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); counts;", {HistType::kTH1D, {axisPt}}); /// Σc0,++
    registry.add("Data/hEtaSoftPiSc0", "#pi^{-} #leftarrow #Sigma_{c}^{0} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                                    /// Σc0
    registry.add("Data/hEtaSoftPiScPlusPlus", "#pi^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                            /// Σc++
    registry.add("Data/hEtaSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0, ++} candidates; #eta ; counts;", {HistType::kTH1D, {axisEta}});                                                      /// Σc0,++
    registry.add("Data/hPhiSoftPiSc0", "#pi^{-} #leftarrow #Sigma_{c}^{0} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                                    /// Σc0
    registry.add("Data/hPhiSoftPiScPlusPlus", "#pi^{+} #leftarrow #Sigma_{c}^{++} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                            /// Σc++
    registry.add("Data/hPhiSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0, ++} candidates; #Phi ; counts;", {HistType::kTH1D, {axisPhi}});                                                      /// Σc0,++

    /// THn for candidate Λc+ and Σc0,++ cut variation
    if (enableTHn) {
      const AxisSpec thnAxisChargeSigmaC{configAxisChargeSigmaC, "charge of SigmaC"};
      const AxisSpec thnAxisMassLambdaC{configAxisMassLambdaC, "inv. mass (K0s p) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPtLambdaC{thnConfigAxisPt, "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtSigmaC{thnConfigAxisPt, "#it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c})"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length #Lambda_{c}^{+} (cm)"};
      const AxisSpec thnAxisDecLengthXY{thnConfigAxisDecLengthXY, "decay length XY #Lambda_{c}^{+} (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle #Lambda_{c}^{+}"};
      const AxisSpec thnAxisCPAXY{thnConfigAxisCPAXY, "cosine of pointing angle XY #Lambda_{c}^{+}"};
      registry.add("hnSigmaC0PlusPlus", "THn for Sigmac from Cascade channel", HistType::kTHnSparseF, {thnAxisChargeSigmaC, thnAxisMassLambdaC, thnAxisPtLambdaC, axisDeltaMassSigmaC, thnAxisDecLength, thnAxisDecLengthXY, thnAxisCPA, thnAxisCPAXY, thnAxisPtSigmaC});
    }
  }

  void processSigmacToLcPi(aod::HfCandScCascades const& candScs,
                           RecoLc const&,
                           aod::Tracks const&)
  {
    for (const auto& candSc : candScs) {
      const auto& candidateLc = candSc.prongLc_as<RecoLc>();
      float ptSc(candSc.pt()), ptLc(candidateLc.pt());
      float const etaSc(candSc.eta()) /*, etaLc(candidateLc.eta())*/;
      float const phiSc(candSc.phi()) /*, phiLc(candidateLc.phi())*/;
      float ptSoftPi(candSc.prong1().pt()), etaSoftPi(candSc.prong1().eta()), phiSoftPi(candSc.prong1().phi());
      double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
      float cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
      float y(-1.);

      auto massLc = HfHelper::invMassLcToK0sP(candidateLc);
      auto massSc = HfHelper::invMassScRecoLcToK0sP(candSc, candidateLc);
      auto deltaMass = massSc - massLc;
      if (candSc.charge() == 0) {
        y = HfHelper::ySc0(candSc);
      } else if (candSc.charge() == 2) {
        y = HfHelper::yScPlusPlus(candSc);
      }
      registry.fill(HIST("Data/hDeltaMassSc0PlusPlus"), deltaMass);           /// Σc(0,++) for both charges
      registry.fill(HIST("Data/hDeltaMassSc0PlusPlusVsPt"), deltaMass, ptSc); /// Σc(0,++) for both charges
      registry.fill(HIST("Data/hEtaSc0PlusPlus"), etaSc);                     /// Σc(0,++) for both charges
      registry.fill(HIST("Data/hEtaSc0PlusPlusVsPt"), etaSc, ptSc);           /// Σc(0,++) for both charges
      registry.fill(HIST("Data/hYSc0PlusPlusVsPt"), y, ptSc);                 /// Σc(0,++) for both charges
      registry.fill(HIST("Data/hPhiSc0PlusPlus"), phiSc);                     /// Σc(0,++) for both charges

      /// fill histograms for softpion from Σc(0,++)
      registry.fill(HIST("Data/hPtSoftPiSc0PlusPlus"), ptSoftPi);
      registry.fill(HIST("Data/hEtaSoftPiSc0PlusPlus"), etaSoftPi);
      registry.fill(HIST("Data/hPhiSoftPiSc0PlusPlus"), phiSoftPi); // π ← Σc0,++

      if (candSc.charge() == 0) {
        registry.fill(HIST("Data/hDeltaMassSc0"), deltaMass);           /// Σc0
        registry.fill(HIST("Data/hDeltaMassSc0VsPt"), deltaMass, ptSc); /// Σc0
        registry.fill(HIST("Data/hEtaSc0"), etaSc);                     /// Σc0
        registry.fill(HIST("Data/hPhiSc0"), phiSc);                     /// Σc0

        /// fill histograms for softpion
        registry.fill(HIST("Data/hPtSoftPiSc0"), ptSoftPi);
        registry.fill(HIST("Data/hEtaSoftPiSc0"), etaSoftPi);
        registry.fill(HIST("Data/hPhiSoftPiSc0"), phiSoftPi); // π ← Σc0

        if (enableTHn) {
          registry.get<THnSparse>(HIST("hnSigmaC0PlusPlus"))->Fill(candSc.charge(), massLc, ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, ptSc);
        }
        /// Σc0
      } else if (candSc.charge() == 2) {
        registry.fill(HIST("Data/hDeltaMassScPlusPlus"), deltaMass);           /// Σc++
        registry.fill(HIST("Data/hDeltaMassScPlusPlusVsPt"), deltaMass, ptSc); /// Σc++
        registry.fill(HIST("Data/hEtaScPlusPlus"), etaSc);                     /// Σc++
        registry.fill(HIST("Data/hPhiScPlusPlus"), phiSc);                     /// Σc++

        /// fill histograms for softpion
        registry.fill(HIST("Data/hPtSoftPiScPlusPlus"), ptSoftPi);
        registry.fill(HIST("Data/hEtaSoftPiScPlusPlus"), etaSoftPi);
        registry.fill(HIST("Data/hPhiSoftPiScPlusPlus"), phiSoftPi); // π ← Σc++

        if (enableTHn) {
          registry.get<THnSparse>(HIST("hnSigmaC0PlusPlus"))->Fill(candSc.charge(), massLc, ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, ptSc);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskSigmacToCascade, processSigmacToLcPi, "Process Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSigmacToCascade>(cfgc)};
}
