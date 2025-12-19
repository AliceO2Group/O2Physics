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

/// \file taskSigmac.cxx
/// \brief Task for Σc0,++ → Λc+(→pK-π+) π-,+ analysis
/// \note Σc0,++ candidates built in O2Physics/PWGHF/TableProducer/HFCandidateCreatorScZeroPlusPlus.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/Utils/utilsSigmac.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TPDGCode.h>

#include <Rtypes.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

struct HfTaskSigmac {
  /// One value of rapidity only
  /// Remember that in Run2 the distinction among GenLimAcc, GenAccMother, GenAcc was done, where:
  ///  - GenLimAcc: Sc in |y|<0.5
  ///  - GenAccMother: Sc in the y-range of the reconstruction ("fiducial acceptance")
  ///  - GenAcc: Sc and Lc in fiducial acceptance, L daughters in acceptance
  /// Properly normalize your results to provide a cross section
  /// OR
  /// consider the new parametrization of the fiducial acceptance (to be seen for reco signal in MC)
  Configurable<float> yCandGenMax{"yCandGenMax", -1, "Maximum generated Sc rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", -1, "Maximum Sc candidate rapidity"};
  Configurable<bool> enableTHn{"enableTHn", false, "enable the usage of THn for Λc+ and Σc0,++"};
  Configurable<bool> addSoftPiDcaToSigmacSparse{"addSoftPiDcaToSigmacSparse", false, "enable the filling of sof-pion dcaXY, dcaZ in the Σc0,++ THnSparse"};
  Configurable<float> deltaMassSigmacRecoMax{"deltaMassSigmacRecoMax", 1000, "Maximum allowed value for Sigmac deltaMass. Conceived to reduce the output size (i.e. reject background above a certain threshold)"};

  bool isMc{};
  static constexpr std::size_t NDaughters{2u};

  using RecoLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;

  /// THn for candidate Λc+ and Σc0,++ cut variation
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {16, 0, 16}, ""};
  ConfigurableAxis thnConfigAxisGenPt{"thnConfigAxisGenPt", {240, 0, 24}, "Gen pt prompt"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {800, 0, 80}, "Gen pt non-prompt"};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisDecLengthXY{"thnConfigAxisDecLengthXY", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisCPAXY{"thnConfigAxisCPAXY", {20, 0.8, 1}, ""};
  ConfigurableAxis configAxisMassLambdaC{"configAxisMassLambdaC", {600, 1.98, 2.58}, ""};
  ConfigurableAxis configAxisDeltaMassSigmaC{"configAxisDeltaMassSigmaC", {200, 0.13, 0.23}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreLcBkg{"thnConfigAxisBdtScoreLcBkg", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreLcNonPrompt{"thnConfigAxisBdtScoreLcNonPrompt", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisSoftPiAbsDca{"thnConfigAxisSoftPiAbsDca", {14, 0., 0.07}, ""};

  /// analysis histograms
  HistogramRegistry registry{
    "registry",
    {/// Σc0
     {"Data/hPtSc0", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSc0", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSc0", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     //{"Data/hDeltaMassSc0", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                             /// Σc++
     {"Data/hPtScPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaScPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiScPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     //{"Data/hDeltaMassScPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}},
     {"Data/hPtSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                                      /// Σc0,++
     {"Data/hPtSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     //{"Data/hDeltaMassSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                                           /// Λc+ ← Σc0
     {"Data/hPtLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     //{"Data/hDeltaMassLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}},
     /// Λc+ ← Σc++
     {"Data/hPtLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     //{"Data/hDeltaMassLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}},
     /// Λc+ ← Σc0,++
     {"Data/hPtLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}}};
  //{"Data/hDeltaMassLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}}}}};

  /// @brief init function, to define the additional analysis histograms
  /// @param
  void init(InitContext&)
  {

    /// To be considered in the future, let's keep the possibility to run in MC also with "data-like" mode (just for TH1 objects)
    // std::array<int, 4> processes {doprocessDataWoMl, doprocessDataWithMl, doprocessMcWoMl, doprocessMcWithMl};
    // if( std::accumulate(processes.begin(), processes.end(), 0) != 1 ) {
    //   LOG(fatal) << "One and only one process function must be enabled. Fix it!";
    // }

    // avoid 2 enabled process functions on data
    if (doprocessDataWoMl && doprocessDataWithMl) {
      LOG(fatal) << "processDataWoMl and processDataWithMl both enabled. Fix it!";
    }
    // avoid 2 enabled process functions on MC
    if (doprocessMcWoMl && doprocessMcWithMl) {
      LOG(fatal) << "processMcWoMl and processMcWithMl both enabled. Fix it!";
    }
    // avoid that in data no ML is used while in MC yes, and viceversa
    if ((doprocessDataWithMl && doprocessMcWoMl) || (doprocessDataWoMl && doprocessMcWithMl)) {
      LOG(fatal) << "process functions with and w/o ML enabled not consistently between data and MC. Fix it! processDataWoMl: " << doprocessDataWoMl << "processDataWithMl: " << doprocessDataWithMl << "processMcWoMl: " << doprocessMcWoMl << "processMcWithMl: " << doprocessMcWithMl;
    }

    /// establish if the analysis is done on Data or MC
    isMc = doprocessMcWoMl || doprocessMcWithMl;

    const AxisSpec thnAxisMassLambdaC{configAxisMassLambdaC, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPtLambdaC{thnConfigAxisPt, "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})"};
    const AxisSpec thnAxisPtSigmaC{thnConfigAxisPt, "#it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c})"};
    const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length #Lambda_{c}^{+} (cm)"};
    const AxisSpec thnAxisDecLengthXY{thnConfigAxisDecLengthXY, "decay length XY #Lambda_{c}^{+} (cm)"};
    const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle #Lambda_{c}^{+}"};
    const AxisSpec thnAxisCPAXY{thnConfigAxisCPAXY, "cosine of pointing angle XY #Lambda_{c}^{+}"};
    const AxisSpec thnAxisOriginMc{3, -0.5, 2.5, "0: none, 1: prompt, 2: non-prompt"};
    const AxisSpec thnAxisChargeSigmaC{3, -0.5, 2.5, "#Sigma_{c}-baryon charge"};
    const AxisSpec thnAxisChannel{4, -0.5, 3.5, "0: direct  1,2,3: resonant"};
    const AxisSpec thnAxisBdtScoreLcBkg{thnConfigAxisBdtScoreLcBkg, "BDT bkg score (Lc)"};
    const AxisSpec thnAxisBdtScoreLcNonPrompt{thnConfigAxisBdtScoreLcNonPrompt, "BDT non-prompt score (Lc)"};
    const AxisSpec thnAxisGenPtLambdaC{thnConfigAxisGenPt, "#it{p}_{T}^{gen}(#Lambda_{c}^{+}) (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtSigmaC{thnConfigAxisGenPt, "#it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtLambdaCBMother{thnConfigAxisGenPtB, "#it{p}_{T}^{gen}(#Lambda_{c}^{+} B mother) (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtSigmaCBMother{thnConfigAxisGenPtB, "#it{p}_{T}^{gen}(#Sigma_{c}^{0,++} B mother) (GeV/#it{c})"};
    const AxisSpec thnAxisSoftPiAbsDcaXY{thnConfigAxisSoftPiAbsDca, "|dca_{xy}|(#pi^{-,+} #leftarrow #Sigma_{c}^{0,++}) (cm)"};
    const AxisSpec thnAxisSoftPiAbsDcaZ{thnConfigAxisSoftPiAbsDca, "|dca_{z}|(#pi^{-,+} #leftarrow #Sigma_{c}^{0,++}) (cm)"};
    const AxisSpec thnAxisGenSigmaCSpecies = {o2::aod::hf_cand_sigmac::Species::NSpecies, -0.5f, +o2::aod::hf_cand_sigmac::Species::NSpecies - 0.5f, "bin 1: #Sigma_{c}(2455), bin 2: #Sigma_{c}(2520)"};
    const AxisSpec thnAxisSigmaCParticleAntiparticle = {o2::aod::hf_cand_sigmac::Conjugated::NConjugated, -0.5f, +o2::aod::hf_cand_sigmac::Conjugated::NConjugated - 0.5f, "bin 1: particle, bin 2: antiparticle"};
    const AxisSpec axisDeltaMassSigmaC{configAxisDeltaMassSigmaC, "#it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2})"};
    registry.add("Data/hDeltaMassSc0", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    registry.add("Data/hDeltaMassScPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    registry.add("Data/hDeltaMassSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    registry.add("Data/hDeltaMassLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    registry.add("Data/hDeltaMassLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    registry.add("Data/hDeltaMassLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {axisDeltaMassSigmaC, {36, 0., 36.}}});
    if (isMc) {
      /////////////////////
      ///   Generated   ///
      /////////////////////
      /// Generated Σc0 signal
      registry.add("MC/generated/hPtGenSc0Sig", "#Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSc0Sig", "#Sigma_{c}^{0} generated signal; #eta^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSc0Sig", "#Sigma_{c}^{0} generated signal; #varphi^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});     /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                      /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      /// Generated Σc++ signal
      registry.add("MC/generated/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #eta^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #varphi^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});     /// soft π
      registry.add("MC/generated/hEtaGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                      /// soft π
      registry.add("MC/generated/hPhiGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      /// Generated Σc0,++ signal
      registry.add("MC/generated/hPtGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});     /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                      /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      /// Generated Λc+ ← Σc0 signal
      registry.add("MC/generated/hPtGenLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc++ signal
      registry.add("MC/generated/hPtGenLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc0,++ signal
      registry.add("MC/generated/hPtGenLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});

      /////////////////////////
      ///   Reconstructed   ///
      /////////////////////////
      /// Reconstructed Σc0 signal
      registry.add("MC/reconstructed/hPtSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #eta(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #varphi(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0SigPrompt", "#Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0SigNonPrompt", "#Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Σc++ signal
      registry.add("MC/reconstructed/hPtScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt", "#Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt", "#Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Σc0,++ signal
      registry.add("MC/reconstructed/hPtSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #eta(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #varphi(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt", "#Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt", "#Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Λc+ ← Σc0 signal
      registry.add("MC/reconstructed/hPtLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc++ signal
      registry.add("MC/reconstructed/hPtLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc0,++ signal
      registry.add("MC/reconstructed/hPtLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {axisDeltaMassSigmaC, {36, 0., 36.}, {4, -0.5, 3.5}}});
    }

    /// THn for candidate Λc+ and Σc0,++ cut variation
    if (enableTHn) {
      std::vector<AxisSpec> axesLambdaCWithMl = {thnAxisPtLambdaC, thnAxisMassLambdaC, thnAxisBdtScoreLcBkg, thnAxisBdtScoreLcNonPrompt, thnAxisOriginMc, thnAxisChannel};
      std::vector<AxisSpec> axesSigmaCWithMl = {thnAxisPtLambdaC, axisDeltaMassSigmaC, thnAxisBdtScoreLcBkg, thnAxisBdtScoreLcNonPrompt, thnAxisOriginMc, thnAxisChannel, thnAxisPtSigmaC, thnAxisChargeSigmaC};
      std::vector<AxisSpec> axesLambdaCWoMl = {thnAxisPtLambdaC, thnAxisMassLambdaC, thnAxisDecLength, thnAxisDecLengthXY, thnAxisCPA, thnAxisCPAXY, thnAxisOriginMc, thnAxisChannel};
      std::vector<AxisSpec> axesSigmaCWoMl = {thnAxisPtLambdaC, axisDeltaMassSigmaC, thnAxisDecLength, thnAxisDecLengthXY, thnAxisCPA, thnAxisCPAXY, thnAxisOriginMc, thnAxisChannel, thnAxisPtSigmaC, thnAxisChargeSigmaC};
      if (isMc) {
        registry.add("MC/generated/hnLambdaCGen", "THn for Lambdac gen", HistType::kTHnSparseF, {thnAxisGenPtLambdaC, thnAxisGenPtLambdaCBMother, thnAxisOriginMc, thnAxisChannel});
        registry.add("MC/generated/hnSigmaCGen", "THn for Sigmac gen", HistType::kTHnSparseF, {thnAxisGenPtSigmaC, thnAxisGenPtSigmaCBMother, thnAxisOriginMc, thnAxisChannel, thnAxisGenPtLambdaC, thnAxisChargeSigmaC, thnAxisGenSigmaCSpecies, thnAxisSigmaCParticleAntiparticle});
        if (doprocessMcWithMl) {
          axesLambdaCWithMl.push_back(thnAxisGenPtLambdaCBMother);
          axesSigmaCWithMl.push_back(thnAxisGenPtSigmaCBMother);
          axesSigmaCWithMl.push_back(thnAxisGenSigmaCSpecies);
          axesSigmaCWithMl.push_back(thnAxisSigmaCParticleAntiparticle);
          if (addSoftPiDcaToSigmacSparse) {
            axesSigmaCWithMl.push_back(thnAxisSoftPiAbsDcaXY);
            axesSigmaCWithMl.push_back(thnAxisSoftPiAbsDcaZ);
          }
          registry.add("hnLambdaC", "THn for Lambdac", HistType::kTHnSparseF, axesLambdaCWithMl);
          registry.add("hnSigmaC", "THn for Sigmac", HistType::kTHnSparseF, axesSigmaCWithMl);
        } else {
          axesLambdaCWoMl.push_back(thnAxisGenPtLambdaCBMother);
          axesSigmaCWoMl.push_back(thnAxisGenPtSigmaCBMother);
          axesSigmaCWoMl.push_back(thnAxisGenSigmaCSpecies);
          axesSigmaCWoMl.push_back(thnAxisSigmaCParticleAntiparticle);
          if (addSoftPiDcaToSigmacSparse) {
            axesSigmaCWoMl.push_back(thnAxisSoftPiAbsDcaXY);
            axesSigmaCWoMl.push_back(thnAxisSoftPiAbsDcaZ);
          }
          registry.add("hnLambdaC", "THn for Lambdac", HistType::kTHnSparseF, axesLambdaCWoMl);
          registry.add("hnSigmaC", "THn for Sigmac", HistType::kTHnSparseF, axesSigmaCWoMl);
        }
      } else {
        if (doprocessDataWithMl) {
          if (addSoftPiDcaToSigmacSparse) {
            axesSigmaCWithMl.push_back(thnAxisSoftPiAbsDcaXY);
            axesSigmaCWithMl.push_back(thnAxisSoftPiAbsDcaZ);
          }
          registry.add("hnLambdaC", "THn for Lambdac", HistType::kTHnSparseF, axesLambdaCWithMl);
          registry.add("hnSigmaC", "THn for Sigmac", HistType::kTHnSparseF, axesSigmaCWithMl);
        } else {
          if (addSoftPiDcaToSigmacSparse) {
            axesSigmaCWoMl.push_back(thnAxisSoftPiAbsDcaXY);
            axesSigmaCWoMl.push_back(thnAxisSoftPiAbsDcaZ);
          }
          registry.add("hnLambdaC", "THn for Lambdac", HistType::kTHnSparseF, axesLambdaCWoMl);
          registry.add("hnSigmaC", "THn for Sigmac", HistType::kTHnSparseF, axesSigmaCWoMl);
        }
      }
    }

  }; /// end init

  /// @brief function to fill the histograms needed in analysis (data)
  /// @param candidatesSc are the reconstructed candidate Σc0,++
  /// @param
  template <bool UseMl, typename CandsLc>
  void fillHistosData(aod::HfCandSc const& candidatesSc,
                      CandsLc const& candidatesLc,
                      aod::Tracks const&)
  {

    /// loop over the candidate Σc0,++
    for (const auto& candSc : candidatesSc) {

      /// rapidity selection on Σc0,++
      /// NB: since in data we cannot tag Sc(2455) and Sc(2520), then we use only Sc(2455) for y selection on reconstructed signal
      if (yCandRecoMax >= 0. && std::abs(HfHelper::ySc0(candSc)) > yCandRecoMax && std::abs(HfHelper::yScPlusPlus(candSc)) > yCandRecoMax) {
        continue;
      }

      const int8_t chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the candidate Σc0,++
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<CandsLc>();
      // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      const int8_t isCandPKPiPiKP = hf_sigmac_utils::isDecayToPKPiToPiKP(candidateLc, candSc);
      double massSc(-1.), massLc(-1.), deltaMass(-1.);
      double ptSc(candSc.pt()), ptLc(candidateLc.pt());
      double etaSc(candSc.eta()), etaLc(candidateLc.eta());
      double phiSc(candSc.phi()), phiLc(candidateLc.phi());
      double ptSoftPi(candSc.prong1().pt()), etaSoftPi(candSc.prong1().eta()), phiSoftPi(candSc.prong1().phi());
      double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
      double cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
      /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
      if (TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PKPi)) {
        massSc = HfHelper::invMassScRecoLcToPKPi(candSc, candidateLc);
        massLc = HfHelper::invMassLcToPKPi(candidateLc);
        deltaMass = massSc - massLc;

        if (deltaMass > deltaMassSigmacRecoMax) {
          /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
          continue;
        }

        /// fill the histograms
        if (chargeSc == o2::aod::hf_cand_sigmac::ChargeNull) {
          registry.fill(HIST("Data/hPtSc0"), ptSc);
          registry.fill(HIST("Data/hEtaSc0"), etaSc);
          registry.fill(HIST("Data/hPhiSc0"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0"), deltaMass, ptSc); // Σc0
          registry.fill(HIST("Data/hPtSoftPiSc0"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtSc0PlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaSc0PlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiSc0PlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0PlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSc0PlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0PlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0PlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLcFromSc0"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0"), deltaMass, ptLc); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLcFromSc0PlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0PlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0PlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0PlusPlus"), deltaMass, ptLc); // Λc+ ← Σc0,++
        } else {                                                                    /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSc0PlusPlus.cxx
          registry.fill(HIST("Data/hPtScPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScPlusPlus"), deltaMass, ptSc); // Σc++
          registry.fill(HIST("Data/hPtSoftPiScPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtSc0PlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaSc0PlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiSc0PlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0PlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSc0PlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0PlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0PlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLcFromScPlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromScPlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromScPlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromScPlusPlus"), deltaMass, ptLc); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLcFromSc0PlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0PlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0PlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0PlusPlus"), deltaMass, ptLc); // Λc+ ← Σc0,++
        }
        /// THn for candidate Σc0,++ cut variation
        if (enableTHn) {
          if (!isMc) {
            /// fill it only if no MC operations are enabled, otherwise fill it in the processMC with the right origin and channel!
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPKPi().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPKPi()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPKPi()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), 0, 0, ptSc, std::abs(chargeSc), softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), 0, 0, ptSc, std::abs(chargeSc));
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0, ptSc, std::abs(chargeSc), softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0, ptSc, std::abs(chargeSc));
              }
            }
          }
        }
      } /// end candidate Λc+ → pK-π+ (and charge conjugate)
      /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
      if (TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PiKP)) {
        massSc = HfHelper::invMassScRecoLcToPiKP(candSc, candidateLc);
        massLc = HfHelper::invMassLcToPiKP(candidateLc);
        deltaMass = massSc - massLc;

        if (deltaMass > deltaMassSigmacRecoMax) {
          /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
          continue;
        }

        /// fill the histograms
        if (chargeSc == o2::aod::hf_cand_sigmac::ChargeNull) {
          registry.fill(HIST("Data/hPtSc0"), ptSc);
          registry.fill(HIST("Data/hEtaSc0"), etaSc);
          registry.fill(HIST("Data/hPhiSc0"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0"), deltaMass, ptSc); // Σc0
          registry.fill(HIST("Data/hPtSoftPiSc0"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtSc0PlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaSc0PlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiSc0PlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0PlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSc0PlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0PlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0PlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLcFromSc0"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0"), deltaMass, ptLc); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLcFromSc0PlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0PlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0PlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0PlusPlus"), deltaMass, ptLc); // Λc+ ← Σc0,++
        } else {                                                                    /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSc0PlusPlus.cxx
          registry.fill(HIST("Data/hPtScPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScPlusPlus"), deltaMass, ptSc); // Σc++
          registry.fill(HIST("Data/hPtSoftPiScPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtSc0PlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaSc0PlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiSc0PlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassSc0PlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSc0PlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSc0PlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSc0PlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLcFromScPlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromScPlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromScPlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromScPlusPlus"), deltaMass, ptLc); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLcFromSc0PlusPlus"), ptLc);
          registry.fill(HIST("Data/hEtaLcFromSc0PlusPlus"), etaLc);
          registry.fill(HIST("Data/hPhiLcFromSc0PlusPlus"), phiLc);
          registry.fill(HIST("Data/hDeltaMassLcFromSc0PlusPlus"), deltaMass, ptLc); // Λc+ ← Σc0,++
        }
        /// THn for candidate Σc0,++ cut variation
        if (enableTHn) {
          if (!isMc) {
            /// fill it only if no MC operations are enabled, otherwise fill it in the processMC with the right origin and channel!
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPiKP().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPiKP()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPiKP()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), 0, 0, ptSc, std::abs(chargeSc), softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), 0, 0, ptSc, std::abs(chargeSc));
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0, ptSc, std::abs(chargeSc), softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0, ptSc, std::abs(chargeSc));
              }
            }
          }
        }
      } /// end candidate Λc+ → π+K-p (and charge conjugate)
    } /// end loop over the candidate Σc0,++

    /// THn for candidate Λc+ cut variation w/o Σc0,++ mass-window cut
    if (enableTHn) {
      /// fill it only if no MC operations are enabled, otherwise fill it in the processMC with the right origin and channel!
      if (!isMc) {
        /// loop over Λc+ candidates w/o Σc0,++ mass-window cut
        for (const auto& candidateLc : candidatesLc) {
          double massLc(-1.);
          double const ptLc(candidateLc.pt());
          double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
          double cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
          if (candidateLc.isSelLcToPKPi() >= 1) {
            massLc = HfHelper::invMassLcToPKPi(candidateLc);
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPKPi().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPKPi()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPKPi()[2]; /// non-prompt score
              }
              registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, outputMl.at(0), outputMl.at(1), 0, 0);
            } else {
              /// fill w/o BDT information
              registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0);
            }
          }
          if (candidateLc.isSelLcToPiKP() >= 1) {
            massLc = HfHelper::invMassLcToPiKP(candidateLc);
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPiKP().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPiKP()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPiKP()[2]; /// non-prompt score
              }
              registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, outputMl.at(0), outputMl.at(1), 0, 0);
            } else {
              /// fill w/o BDT information
              registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, 0, 0);
            }
          }
        }
      }
    } /// end THn for candidate Λc+ cut variation w/o Σc0,++ mass-window cut
  }; /// end fillHistosData

  /// @brief function to remap the value of the resonant decay channel to fit the binning of the thnAxisChannel axis
  /// @param channel the value obtained from candidateLc.flagMcDecayChanGen() or particleLc.flagMcDecayChanGen()
  int remapResoChannelLc(int channel)
  {
    switch (channel) {
      case 0:
        // direct channel
        return 0;
      case o2::hf_decay::hf_cand_3prong::DecayChannelResonant::LcToPKstar0:
        return 1;
      case o2::hf_decay::hf_cand_3prong::DecayChannelResonant::LcToDeltaplusplusK:
        return 2;
      case o2::hf_decay::hf_cand_3prong::DecayChannelResonant::LcToL1520Pi:
        return 3;
    }
    return -1;
  }

  /// @brief function to fill the histograms needed in analysis (MC)
  /// @param candidatesSc are the reconstructed candidate Σc0,++ with MC info
  /// @param mcParticles are the generated particles with flags wheter they are Σc0,++ or not
  /// @param
  template <bool UseMl, typename CandsLc>
  void fillHistosMc(soa::Join<aod::HfCandSc, aod::HfCandScMcRec> const& candidatesSc,
                    soa::Join<aod::McParticles, aod::HfCandScMcGen> const& mcParticlesSc,
                    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticlesLc,
                    aod::McParticles const& mcParticles, // this establishes the type of particle obtained with the .mcParticle() getter
                    CandsLc const& candidatesLc,
                    aod::TracksWMc const&)
  {

    /// loop over Sc generated particles
    for (const auto& particle : mcParticlesSc) {

      /// reject immediately particles different from Σc0,++
      bool const isSc0Gen = (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::Sc0ToPKPiPi);
      bool const isScStar0Gen = (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStar0ToPKPiPi);
      bool const isScPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScplusplusToPKPiPi);
      bool const isScStarPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStarPlusPlusToPKPiPi);
      if (!isSc0Gen && !isScPlusPlusGen && !isScStar0Gen && !isScStarPlusPlusGen) {
        continue;
      }

      /// look for generated particles in acceptance
      /*
         One value of rapidity only
         Remember that in Run2 the distinction among GenLimAcc, GenAccMother, GenAcc was done, where:
          - GenLimAcc: Sc in |y|<0.5
          - GenAccMother: Sc in the y-range of the reconstruction ("fiducial acceptance")
          - GenAcc: Sc and Lc in fiducial acceptance, L daughters in acceptance
         Properly normalize your results to provide a cross section
         OR
         consider the new parametrization of the fiducial acceptance (to be seen for reco signal in MC)
      */
      if (yCandGenMax >= 0.) {
        double mass = -1;
        if (isSc0Gen) {
          mass = o2::constants::physics::MassSigmaC0;
        } else if (isScPlusPlusGen) {
          mass = o2::constants::physics::MassSigmaCPlusPlus;
        } else if (isScStar0Gen) {
          mass = o2::constants::physics::MassSigmaCStar0;
        } else if (isScStarPlusPlusGen) {
          mass = o2::constants::physics::MassSigmaCStarPlusPlus;
        }
        if (mass > -1. && std::abs(RecoDecay::y(particle.pVector(), mass)) > yCandGenMax) {
          continue;
        }
      }

      /// Get the kinematic information of Σc0,++ and the daughters
      /// Get information about origin (prompt, non-prompt)
      /// Get information about decay Λc+ channel (direct, resonant)
      double ptGenSc(particle.pt()), etaGenSc(particle.eta()), phiGenSc(particle.phi());
      double ptGenScBMother(-1.);
      auto arrayDaughtersIds = particle.daughtersIds();
      if (arrayDaughtersIds.size() != NDaughters) {
        /// This should never happen
        LOG(fatal) << "generated Σc0,++ has a number of daughter particles different than 2";
        continue;
      }
      double ptGenLc(-1.), ptGenSoftPi(-1.);
      double etaGenLc(-1.), etaGenSoftPi(-1.);
      double phiGenLc(-1.), phiGenSoftPi(-1.);
      int origin = -1;
      int8_t channel = -1;
      auto daughter0 = mcParticles.rawIteratorAt(arrayDaughtersIds[0]);
      if (std::abs(daughter0.pdgCode()) == o2::constants::physics::Pdg::kLambdaCPlus) {
        /// daughter 0 is the Λc+, daughter 1 the soft π
        auto daugLc = mcParticlesLc.rawIteratorAt(arrayDaughtersIds[0]);
        auto daugSoftPi = mcParticles.rawIteratorAt(arrayDaughtersIds[1]);
        ptGenLc = daugLc.pt();
        etaGenLc = daugLc.eta();
        phiGenLc = daugLc.phi();
        origin = daugLc.originMcGen();
        channel = daugLc.flagMcDecayChanGen();
        ptGenSoftPi = daugSoftPi.pt();
        etaGenSoftPi = daugSoftPi.eta();
        phiGenSoftPi = daugSoftPi.phi();
      } else if (std::abs(daughter0.pdgCode()) == kPiPlus) {
        /// daughter 0 is the soft π, daughter 1 the Λc+
        auto daugLc = mcParticlesLc.rawIteratorAt(arrayDaughtersIds[1]);
        auto daugSoftPi = mcParticles.rawIteratorAt(arrayDaughtersIds[0]);
        ptGenLc = daugLc.pt();
        etaGenLc = daugLc.eta();
        phiGenLc = daugLc.phi();
        origin = daugLc.originMcGen();
        channel = daugLc.flagMcDecayChanGen();
        ptGenSoftPi = daugSoftPi.pt();
        etaGenSoftPi = daugSoftPi.eta();
        phiGenSoftPi = daugSoftPi.phi();
      }
      channel = remapResoChannelLc(channel);

      /// Fill histograms
      int sigmacSpecies = -1;
      if (isSc0Gen || isScPlusPlusGen) {
        sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2455;
      } else if (isScStar0Gen || isScStarPlusPlusGen) {
        sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2520;
      }
      if (isSc0Gen || isScStar0Gen) {
        /// Generated Σc0 and Λc+ ← Σc0 signals
        registry.fill(HIST("MC/generated/hPtGenSc0Sig"), ptGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSc0Sig"), etaGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSc0Sig"), phiGenSc, origin, channel); /// Generated Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSc0Sig"), ptGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSc0Sig"), etaGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSc0Sig"), phiGenSoftPi, origin, channel); /// Generated π ← Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenLcFromSc0Sig"), ptGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenLcFromSc0Sig"), etaGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenLcFromSc0Sig"), phiGenLc, origin, channel); /// Generated Λc+ ← Σc0 signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSc0PlusPlusSig"), etaGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSc0PlusPlusSig"), phiGenSc, origin, channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSc0PlusPlusSig"), etaGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSc0PlusPlusSig"), phiGenSoftPi, origin, channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenLcFromSc0PlusPlusSig"), etaGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenLcFromSc0PlusPlusSig"), phiGenLc, origin, channel); /// Generated Λc+ ← Σc0,++ signal
        int8_t const particleAntiparticle = particle.particleAntiparticle();
        if (origin == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/hnSigmaCGen"), ptGenSc, ptGenScBMother, origin, channel, ptGenLc, 0, sigmacSpecies, particleAntiparticle);
        } else {
          ptGenScBMother = mcParticlesSc.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("MC/generated/hnSigmaCGen"), ptGenSc, ptGenScBMother, origin, channel, ptGenLc, 0, sigmacSpecies, particleAntiparticle);
        }

        // debug -- uncomment if needed
        // it should be solved after the implementation of ev. selection for generated SigmaC particles
        // if(origin != RecoDecay::OriginType::Prompt && origin != RecoDecay::OriginType::NonPrompt) {
        //  LOG(info) << "   --> (Sc0 gen) origin " << static_cast<int>(origin) << ", particle.originMcGen() " << static_cast<int>(particle.originMcGen()) << ", particle.flagMcMatchGen() " << static_cast<int>(particle.flagMcMatchGen()) << ", pdg " << particle.pdgCode();
        //}

      } else if (isScPlusPlusGen || isScStarPlusPlusGen) {
        /// Generated Σc++ and Λc+ ← Σc++ signals
        registry.fill(HIST("MC/generated/hPtGenScPlusPlusSig"), ptGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenScPlusPlusSig"), etaGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenScPlusPlusSig"), phiGenSc, origin, channel); /// Generated Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiScPlusPlusSig"), etaGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiScPlusPlusSig"), phiGenSoftPi, origin, channel); /// Generated π ← Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenLcFromScPlusPlusSig"), ptGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenLcFromScPlusPlusSig"), etaGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenLcFromScPlusPlusSig"), phiGenLc, origin, channel); /// Generated Λc+ ← Σc++ signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSc0PlusPlusSig"), etaGenSc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSc0PlusPlusSig"), phiGenSc, origin, channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSc0PlusPlusSig"), etaGenSoftPi, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSc0PlusPlusSig"), phiGenSoftPi, origin, channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hEtaGenLcFromSc0PlusPlusSig"), etaGenLc, origin, channel);
        registry.fill(HIST("MC/generated/hPhiGenLcFromSc0PlusPlusSig"), phiGenLc, origin, channel); /// Generated Λc+ ← Σc0,++ signal
        int8_t const particleAntiparticle = particle.particleAntiparticle();
        if (origin == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/hnSigmaCGen"), ptGenSc, ptGenScBMother, origin, channel, ptGenLc, 2, sigmacSpecies, particleAntiparticle);
        } else {
          ptGenScBMother = mcParticlesSc.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("MC/generated/hnSigmaCGen"), ptGenSc, ptGenScBMother, origin, channel, ptGenLc, 2, sigmacSpecies, particleAntiparticle);
        }

        // debug -- uncomment if needed
        // it should be solved after the implementation of ev. selection for generated SigmaC particles
        // if(origin != RecoDecay::OriginType::Prompt && origin != RecoDecay::OriginType::NonPrompt) {
        //  LOG(info) << "   --> (Sc++ gen) origin " << static_cast<int>(origin) << ", particle.originMcGen() " << static_cast<int>(particle.originMcGen()) << ", particle.flagMcMatchGen() " << static_cast<int>(particle.flagMcMatchGen()) << ", pdg " << particle.pdgCode();
        //}
      }

    } /// end loop over Sc generated particles

    /// loop over Lc generated particles
    for (const auto& particle : mcParticlesLc) {
      if (std::abs(particle.flagMcMatchGen()) != hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
        continue;
      }
      if (yCandGenMax >= 0. && std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus)) > yCandGenMax) {
        continue;
      }
      double ptGenLc(particle.pt()), ptGenLcBMother(-1.);
      int const origin = particle.originMcGen();
      int channel = particle.flagMcDecayChanGen();
      channel = remapResoChannelLc(channel);
      if (origin == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("MC/generated/hnLambdaCGen"), ptGenLc, ptGenLcBMother, origin, channel);
      } else {
        ptGenLcBMother = mcParticlesLc.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("MC/generated/hnLambdaCGen"), ptGenLc, ptGenLcBMother, origin, channel);
      }
    } /// end loop over Lc generated particles

    /// reconstructed Σc0,++ matched to MC
    for (const auto& candSc : candidatesSc) {

      /// rapidity selection on Σc0,++
      /// NB: since in data we cannot tag Sc(2455) and Sc(2520), then we use only Sc(2455) for y selection on reconstructed signal
      if (yCandRecoMax >= 0. && std::abs(HfHelper::ySc0(candSc)) > yCandRecoMax && std::abs(HfHelper::yScPlusPlus(candSc)) > yCandRecoMax) {
        continue;
      }

      /// electric charge
      const int8_t chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the Σc0
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<CandsLc>();
      const int8_t isCandPKPiPiKP = hf_sigmac_utils::isDecayToPKPiToPiKP(candidateLc, candSc);

      // candidateLc.flagMcDecayChanRec();

      bool const isTrueSc0Reco = std::abs(candSc.flagMcMatchRec()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::Sc0ToPKPiPi;
      bool const isTrueScStar0Reco = std::abs(candSc.flagMcMatchRec()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStar0ToPKPiPi;
      bool const isTrueScPlusPlusReco = std::abs(candSc.flagMcMatchRec()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScplusplusToPKPiPi;
      bool const isTrueScStarPlusPlusReco = std::abs(candSc.flagMcMatchRec()) == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStarPlusPlusToPKPiPi;
      if (!isTrueSc0Reco && !isTrueScStar0Reco && !isTrueScPlusPlusReco && !isTrueScStarPlusPlusReco) {
        continue;
      }
      int sigmacSpecies = -1;

      /// debug
      if ((isTrueSc0Reco || isTrueScStar0Reco) && chargeSc != o2::aod::hf_cand_sigmac::ChargeNull) {
        /// this should never happen
        LOG(fatal) << "isTrueSc0Reco=" << isTrueSc0Reco << ", isTrueScStar0Reco=" << isTrueScStar0Reco << ", but chargeSc = " << static_cast<int>(chargeSc) << "! Not possible, abort...";
      }
      if ((isTrueScPlusPlusReco || isTrueScStarPlusPlusReco) && std::abs(chargeSc) != o2::aod::hf_cand_sigmac::ChargePlusPlus) {
        /// this should never happen
        LOG(fatal) << "isTrueScPlusPlusReco=" << isTrueScPlusPlusReco << ", isTrueScStarPlusPlusReco=" << isTrueScStarPlusPlusReco << ", but chargeSc = " << static_cast<int>(chargeSc) << "! Not possible, abort...";
      }

      if ((isTrueSc0Reco || isTrueScStar0Reco) && (chargeSc == o2::aod::hf_cand_sigmac::ChargeNull)) {
        /// Reconstructed Σc0 signal
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        int indexMcScRec = -1;
        if (isTrueSc0Reco) {
          // Σc0(2455)
          indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), o2::constants::physics::Pdg::kSigmaC0, true);
          sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2455;
        } else if (isTrueScStar0Reco) {
          // Σc0(2520)
          indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), o2::constants::physics::Pdg::kSigmaCStar0, true);
          sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2520;
        }
        auto particleSc = mcParticles.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = std::array{candidateLc.template prong0_as<aod::TracksWMc>(), candidateLc.template prong1_as<aod::TracksWMc>(), candidateLc.template prong2_as<aod::TracksWMc>()};
        int8_t sign = 0;
        int const indexMcLcRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughtersLc, o2::constants::physics::Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLc = mcParticles.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSc.prong1_as<aod::TracksWMc>().mcParticle();

        // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSc(-1.), massLc(-1.), deltaMass(-1.);
        double ptSc(candSc.pt()), ptLc(candidateLc.pt());
        double etaSc(candSc.eta()), etaLc(candidateLc.eta());
        double phiSc(candSc.phi()), phiLc(candidateLc.phi());
        double ptSoftPi(candSc.prong1_as<aod::TracksWMc>().pt()), etaSoftPi(candSc.prong1_as<aod::TracksWMc>().eta()), phiSoftPi(candSc.prong1_as<aod::TracksWMc>().phi());
        double ptGenSc(particleSc.pt()), ptGenLc(particleLc.pt()), ptGenSoftPi(particleSoftPi.pt());
        double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
        double cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
        int const origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
        channel = remapResoChannelLc(channel);

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PKPi)) && std::abs(candidateLc.template prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kProton) {
          massSc = HfHelper::invMassScRecoLcToPKPi(candSc, candidateLc);
          massLc = HfHelper::invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLc;

          if (deltaMass > deltaMassSigmacRecoMax) {
            /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
            continue;
          }

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSc0Sig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0Sig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0Sig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0Sig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0Sig"), deltaMass, ptSc, channel);
          ////////////////////////////////////////////////////////////////////////////////////////////////////// Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0Sig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0Sig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0Sig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0Sig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0Sig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0Sig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0Sig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0Sig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0Sig"), deltaMass, ptLc, channel);
          ////////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          ////////////////////////////////////////////////////////////////////////////////////////////////////// Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0,++ signal

          if (origin == RecoDecay::OriginType::Prompt) {
            /// prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigPrompt"), deltaMass, ptSc, channel);               // Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt"), deltaMass, ptLc, channel);         // Λc+ ← Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            /// non-prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigNonPrompt"), deltaMass, ptSc, channel);               // Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt"), deltaMass, ptLc, channel);         // Λc+ ← Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          }

          /// THn for candidate Σc0,++ cut variation
          if (enableTHn) {
            int8_t const particleAntiparticle = candSc.particleAntiparticle();
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPKPi().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPKPi()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPKPi()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            }
          }

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PiKP)) && std::abs(candidateLc.template prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = HfHelper::invMassScRecoLcToPiKP(candSc, candidateLc);
          massLc = HfHelper::invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLc;

          if (deltaMass > deltaMassSigmacRecoMax) {
            /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
            continue;
          }

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSc0Sig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0Sig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0Sig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0Sig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0Sig"), deltaMass, ptSc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////////////// Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0Sig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0Sig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0Sig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0Sig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0Sig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0Sig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0Sig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0Sig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0Sig"), deltaMass, ptLc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////////////// Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0,++ signal

          if (origin == RecoDecay::OriginType::Prompt) {
            /// prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigPrompt"), deltaMass, ptSc, channel);               // Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt"), deltaMass, ptLc, channel);         // Λc+ ← Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            /// non-prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigNonPrompt"), deltaMass, ptSc, channel);               // Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt"), deltaMass, ptLc, channel);         // Λc+ ← Σc0 signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          }

          /// THn for candidate Σc0,++ cut variation
          if (enableTHn) {
            int8_t const particleAntiparticle = candSc.particleAntiparticle();
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPiKP().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPiKP()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPiKP()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            }
          }

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
        /// end reconstructed Σc0 signal
      } else if ((isTrueScPlusPlusReco || isTrueScStarPlusPlusReco) && (std::abs(chargeSc) == o2::aod::hf_cand_sigmac::ChargePlusPlus)) {
        /// Reconstructed Σc++ signal
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        int indexMcScRec = -1;
        if (isTrueScPlusPlusReco) {
          // Σc0(2455)
          indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), o2::constants::physics::Pdg::kSigmaCPlusPlus, true);
          sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2455;
        } else if (isTrueScStarPlusPlusReco) {
          // Σc0(2520)
          indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), o2::constants::physics::Pdg::kSigmaCStarPlusPlus, true);
          sigmacSpecies = o2::aod::hf_cand_sigmac::Sc2520;
        }
        auto particleSc = mcParticles.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = std::array{candidateLc.template prong0_as<aod::TracksWMc>(), candidateLc.template prong1_as<aod::TracksWMc>(), candidateLc.template prong2_as<aod::TracksWMc>()};
        int8_t sign = 0;
        int const indexMcLcRec = RecoDecay::getMatchedMCRec<false, false, false, true, true>(mcParticles, arrayDaughtersLc, o2::constants::physics::Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLc = mcParticles.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSc.prong1_as<aod::TracksWMc>().mcParticle();

        // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSc(-1.), massLc(-1.), deltaMass(-1.);
        double ptSc(candSc.pt()), ptLc(candidateLc.pt());
        double etaSc(candSc.eta()), etaLc(candidateLc.eta());
        double phiSc(candSc.phi()), phiLc(candidateLc.phi());
        double ptSoftPi(candSc.prong1_as<aod::TracksWMc>().pt()), etaSoftPi(candSc.prong1_as<aod::TracksWMc>().eta()), phiSoftPi(candSc.prong1_as<aod::TracksWMc>().phi());
        double ptGenSc(particleSc.pt()), ptGenLc(particleLc.pt()), ptGenSoftPi(particleSoftPi.pt());
        double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
        double cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
        int const origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±; FIXME: DecayChannelResonant
        channel = remapResoChannelLc(channel);

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PKPi)) && std::abs(candidateLc.template prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kProton) {
          massSc = HfHelper::invMassScRecoLcToPKPi(candSc, candidateLc);
          massLc = HfHelper::invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLc;

          if (deltaMass > deltaMassSigmacRecoMax) {
            /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
            continue;
          }

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, channel);
          //////////////////////////////////////////////////////////////////////////////////////////////// Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromScPlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromScPlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromScPlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromScPlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig"), deltaMass, ptLc, channel);
          //////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc++ signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          /////////////////////////////////////////////////////////////////////////////////////////////// Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          ////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0,++ signal

          if (origin == RecoDecay::OriginType::Prompt) {
            /// prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, channel);        // Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt"), deltaMass, ptLc, channel);  // Λc+ ← Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            /// non-prompt signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, channel);        // Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt"), deltaMass, ptLc, channel);  // Λc+ ← Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          }

          /// THn for candidate Σc0,++ cut variation
          if (enableTHn) {
            int8_t const particleAntiparticle = candSc.particleAntiparticle();
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPKPi().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPKPi()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPKPi()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            }
          }

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PiKP)) && std::abs(candidateLc.template prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = HfHelper::invMassScRecoLcToPiKP(candSc, candidateLc);
          massLc = HfHelper::invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLc;

          if (deltaMass > deltaMassSigmacRecoMax) {
            /// the reconstructed deltaMass is too large, let's ignore this candidate for TH1 / THnSparse filling
            continue;
          }

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, channel);
          ////////////////////////////////////////////////////////////////////////////////////////////////////// Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromScPlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromScPlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromScPlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromScPlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig"), deltaMass, ptLc, channel);
          ///////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc++ signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          ///////////////////////////////////////////////////////////////////////////////////////////////////// Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          //////////////////////////////////////////////////////////////////////////////////////////////////// Λc+ ← Σc0,++ signal

          if (origin == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, channel);        // Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt"), deltaMass, ptLc, channel);  // Λc+ ← Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, channel);        // Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt"), deltaMass, ptLc, channel);  // Λc+ ← Σc++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel);       // Σc0,++ signal
            registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal
          }

          /// THn for candidate Σc0,++ cut variation
          if (enableTHn) {
            int8_t const particleAntiparticle = candSc.particleAntiparticle();
            const float softPiAbsDcaXY = std::abs(candSc.softPiDcaXY());
            const float softPiAbsDcaZ = std::abs(candSc.softPiDcaZ());
            if constexpr (UseMl) {
              /// fill with ML information
              /// BDT index 0: bkg score; BDT index 2: non-prompt score
              std::array<float, 2> outputMl{-1., -1.};
              if (candidateLc.mlProbLcToPiKP().size() > 0) {
                outputMl.at(0) = candidateLc.mlProbLcToPiKP()[0]; /// bkg score
                outputMl.at(1) = candidateLc.mlProbLcToPiKP()[2]; /// non-prompt score
              }
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, outputMl.at(0), outputMl.at(1), origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            } else {
              /// fill w/o BDT information
              if (addSoftPiDcaToSigmacSparse) {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle, softPiAbsDcaXY, softPiAbsDcaZ);
              } else {
                registry.get<THnSparse>(HIST("hnSigmaC"))->Fill(ptLc, deltaMass, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, ptSc, std::abs(chargeSc), candSc.ptBhadMotherPart(), sigmacSpecies, particleAntiparticle);
              }
            }
          }

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
      } /// end reconstructed Σc++ signal

    } /// end loop on reconstructed Σc0,++

    /// THn for candidate Λc+ cut variation w/o Σc0,++ mass-window cut
    if (enableTHn) {
      /// loop over Λc+ candidates w/o Σc0,++ mass-window cut
      for (const auto& candidateLc : candidatesLc) {
        if (std::abs(candidateLc.flagMcMatchRec()) != hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
          continue;
        }
        double massLc(-1.);
        double const ptLc(candidateLc.pt());
        double decLengthLc(candidateLc.decayLength()), decLengthXYLc(candidateLc.decayLengthXY());
        double cpaLc(candidateLc.cpa()), cpaXYLc(candidateLc.cpaXY());
        int const origin = candidateLc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
        channel = remapResoChannelLc(channel);
        int pdgAbs = -1;
        if (candidateLc.template prong0_as<aod::TracksWMc>().has_mcParticle()) {
          pdgAbs = std::abs(candidateLc.template prong0_as<aod::TracksWMc>().mcParticle().pdgCode());
        }
        if (candidateLc.isSelLcToPKPi() >= 1 && pdgAbs == kProton) {
          massLc = HfHelper::invMassLcToPKPi(candidateLc);
          if constexpr (UseMl) {
            /// fill with ML information
            /// BDT index 0: bkg score; BDT index 2: non-prompt score
            std::array<float, 2> outputMl{-1., -1.};
            if (candidateLc.mlProbLcToPKPi().size() > 0) {
              outputMl.at(0) = candidateLc.mlProbLcToPKPi()[0]; /// bkg score
              outputMl.at(1) = candidateLc.mlProbLcToPKPi()[2]; /// non-prompt score
            }
            registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, outputMl.at(0), outputMl.at(1), origin, channel, candidateLc.ptBhadMotherPart());
          } else {
            /// fill w/o BDT information
            registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, candidateLc.ptBhadMotherPart());
          }
        }
        if (candidateLc.isSelLcToPiKP() >= 1 && pdgAbs == kPiPlus) {
          massLc = HfHelper::invMassLcToPiKP(candidateLc);
          if constexpr (UseMl) {
            /// fill with ML information
            /// BDT index 0: bkg score; BDT index 2: non-prompt score
            std::array<float, 2> outputMl{-1., -1.};
            if (candidateLc.mlProbLcToPiKP().size() > 0) {
              outputMl.at(0) = candidateLc.mlProbLcToPiKP()[0]; /// bkg score
              outputMl.at(1) = candidateLc.mlProbLcToPiKP()[2]; /// non-prompt score
            }
            registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, outputMl.at(0), outputMl.at(1), origin, channel, candidateLc.ptBhadMotherPart());
          } else {
            /// fill w/o BDT information
            registry.get<THnSparse>(HIST("hnLambdaC"))->Fill(ptLc, massLc, decLengthLc, decLengthXYLc, cpaLc, cpaXYLc, origin, channel, candidateLc.ptBhadMotherPart());
          }
        }
      }

    } /// end THn for candidate Λc+ cut variation w/o Σc0,++ mass-window cut

  }; /// end fillHistosMc

  /// @brief process function to fill the histograms needed in analysis w/o ML information (data)
  void processDataWoMl(aod::HfCandSc const& candidatesSc,
                       RecoLc const& candidatesLc,
                       aod::Tracks const& tracks)
  {
    fillHistosData<false>(candidatesSc, candidatesLc, tracks);
  }
  PROCESS_SWITCH(HfTaskSigmac, processDataWoMl, "Process data w/o ML information on Lc", true);

  /// @brief process function to fill the histograms needed in analysis with ML information (data)
  void processDataWithMl(aod::HfCandSc const& candidatesSc,
                         soa::Join<RecoLc, aod::HfMlLcToPKPi> const& candidatesLc,
                         aod::Tracks const& tracks)
  {
    fillHistosData<true>(candidatesSc, candidatesLc, tracks);
  }
  PROCESS_SWITCH(HfTaskSigmac, processDataWithMl, "Process data with ML information on Lc", false);

  /// @brief process function to fill the histograms needed in analysis w/o ML information (MC)
  void processMcWoMl(soa::Join<aod::HfCandSc, aod::HfCandScMcRec> const& candidatesSc,
                     soa::Join<aod::McParticles, aod::HfCandScMcGen> const& mcParticlesSc,
                     soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticlesLc,
                     aod::McParticles const& mcParticles, // this establishes the type of particle obtained with the .mcParticle() getter
                     soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const& candidatesLc,
                     aod::TracksWMc const& tracksWithMc)
  {
    fillHistosMc<false>(candidatesSc, mcParticlesSc, mcParticlesLc, mcParticles, candidatesLc, tracksWithMc);
  }
  PROCESS_SWITCH(HfTaskSigmac, processMcWoMl, "Process MC w/o ML information on Lc", false);

  /// @brief process function to fill the histograms needed in analysis with ML information (MC)
  void processMcWithMl(soa::Join<aod::HfCandSc, aod::HfCandScMcRec> const& candidatesSc,
                       soa::Join<aod::McParticles, aod::HfCandScMcGen> const& mcParticlesSc,
                       soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticlesLc,
                       aod::McParticles const& mcParticles, // this establishes the type of particle obtained with the .mcParticle() getter
                       soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec, aod::HfMlLcToPKPi> const& candidatesLc,
                       aod::TracksWMc const& tracksWithMc)
  {
    fillHistosMc<true>(candidatesSc, mcParticlesSc, mcParticlesLc, mcParticles, candidatesLc, tracksWithMc);
  }
  PROCESS_SWITCH(HfTaskSigmac, processMcWithMl, "Process MC with ML information on Lc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSigmac>(cfgc)};
}
