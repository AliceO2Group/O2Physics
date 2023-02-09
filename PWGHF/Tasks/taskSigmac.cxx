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
/// \note Σc0,++ candidates built in O2Physics/PWGHF/TableProducer/HFCandidateCreatorSigmacZeroPlusPlus.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_cand_sc;

struct HfTaskSigmac {

  /// One value of rapidity only
  /// Remember that in Run2 the distinction among GenLimAcc, GenAccMother, GenAcc was done, where:
  ///  - GenLimAcc: Sc in |y|<0.5
  ///  - GenAccMother: Sc in the y-range of the reconstruction ("fiducial acceptance")
  ///  - GenAcc: Sc and Lc in fiducial acceptance, L daughters in acceptance
  /// Properly normalize your results to provide a cross section
  /// OR
  /// consider the new parametrization of the fiducial acceptance (to be seen for reco signal in MC)
  Configurable<float> yCandMax{"yCandMax", -1, "Sigmac rapidity"};

  /// analysis histograms
  HistogramRegistry registry{
    "registry",
    {/// Σc0
     {"Data/hPtSigmacZero", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSigmacZero", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSigmacZero", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassSigmacZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSigmacZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiSigmacZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiSigmacZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                      /// Σc++
     {"Data/hPtSigmacPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSigmacPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSigmacPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassSigmacPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSigmacPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiSigmacPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiSigmacPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                            /// Σc0,++
     {"Data/hPtSigmacZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSigmacZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSigmacZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassSigmacZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSigmacZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiSigmacZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiSigmacZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                                    /// Λc+ ← Σc0
     {"Data/hPtLambdaCFromSigmacZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromSigmacZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromSigmacZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromSigmacZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc++
     {"Data/hPtLambdaCFromSigmacPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromSigmacPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromSigmacPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromSigmacPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc0,++
     {"Data/hPtLambdaCFromSigmacZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromSigmacZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromSigmacZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromSigmacZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}};

  /// @brief init function, to define the additional analysis histograms
  /// @param
  void init(InitContext&)
  {
    if (doprocessMc) {
      /////////////////////
      ///   Generated   ///
      /////////////////////
      /// Generated Σc0 signal
      registry.add("MC/generated/hPtGenSigmacZeroSig", "#Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSigmacZeroSig", "#Sigma_{c}^{0} generated signal; #eta^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSigmacZeroSig", "#Sigma_{c}^{0} generated signal; #varphi^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Σc++ signal
      registry.add("MC/generated/hPtGenSigmacPlusPlusSig", "#Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSigmacPlusPlusSig", "#Sigma_{c}^{++} generated signal; #eta^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSigmacPlusPlusSig", "#Sigma_{c}^{++} generated signal; #varphi^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Σc0,++ signal
      registry.add("MC/generated/hPtGenSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Λc+ ← Σc0 signal
      registry.add("MC/generated/hPtGenLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc++ signal
      registry.add("MC/generated/hPtGenLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc0,++ signal
      registry.add("MC/generated/hPtGenLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});

      /////////////////////////
      ///   Reconstructed   ///
      /////////////////////////
      /// Reconstructed Σc0 signal
      registry.add("MC/reconstructed/hPtSigmacZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSigmacZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSigmacZeroSig", "#Sigma_{c}^{0} reconstructed signal; #eta(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSigmacZeroSig", "#Sigma_{c}^{0} reconstructed signal; #varphi(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroSigPrompt", "#Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroSigNonPrompt", "#Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSigmacZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Σc++ signal
      registry.add("MC/reconstructed/hPtSigmacPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSigmacPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSigmacPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSigmacPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacPlusPlusSigPrompt", "#Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacPlusPlusSigNonPrompt", "#Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSigmacPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Σc0,++ signal
      registry.add("MC/reconstructed/hPtSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #eta(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #varphi(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigPrompt", "#Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigNonPrompt", "#Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSigmacZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Λc+ ← Σc0 signal
      registry.add("MC/reconstructed/hPtLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc++ signal
      registry.add("MC/reconstructed/hPtLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc0,++ signal
      registry.add("MC/reconstructed/hPtLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
    }
  }; /// end init

  /// @brief Function to determine if the reconstructed candidate Σc0,++ decays into Λc+ → pK-π+, Λc+ → π+K-p or both
  /// @tparam L template for LambdaC daughter of Sigmac candidate
  /// @tparam S template for Sigmac candidate
  /// @param candLambdaC LambdaC daughter of Sigmac candidate
  /// @param candSigmac Sigmac candidate
  /// @return 0: none; 1: only Λc+ → pK-π+ possible; 2: Λc+ → π+K-p possible; 3: both possible
  template <typename L, typename S>
  int isDecayToPKPiToPiKP(L& candLambdaC, S& candSigmac)
  {
    int channel = 0;
    if ((candLambdaC.isSelLcToPKPi() >= 1) && candSigmac.statusSpreadLcMinvPKPiFromPDG()) {
      // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      channel += 1;
    }
    if ((candLambdaC.isSelLcToPiKP() >= 1) && candSigmac.statusSpreadLcMinvPiKPFromPDG()) {
      // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      channel += 2;
    }
    return channel; /// 0: none; 1: pK-π+ only; 2: π+K-p only; 3: both possible
  }

  /// @brief process function to fill the histograms needed in analysis (data)
  /// @param candidatesSigmac are the reconstructed candidate Σc0,++
  /// @param
  void process(const aod::HfCandSigmac& candidatesSigmac,
               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&, const soa::Join<aod::Tracks, aod::TracksDCA>&)
  {

    /// loop over the candidate Σc0,++
    for (auto& candSigmac : candidatesSigmac) {

      const int chargeSigmac = candSigmac.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the candidate Σc0,++
      /// and understand which mass hypotheses are possible
      const auto& candLambdaC = candSigmac.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
      // const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmac.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      // const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmac.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candLambdaC, candSigmac);
      double massSigmac(-1.), massLambdaC(-1.), deltaMass(-1.);
      double ptSigmac(candSigmac.pt()), ptLambdaC(candLambdaC.pt());
      double etaSigmac(candSigmac.eta()), etaLambdaC(candLambdaC.eta());
      double phiSigmac(candSigmac.phi()), phiLambdaC(candLambdaC.phi());
      double ptSoftPi(candSigmac.prong1().pt()), etaSoftPi(candSigmac.prong1().eta()), phiSoftPi(candSigmac.prong1().phi());
      /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
        massSigmac = invMassScRecoLcToPKPi(candSigmac);
        massLambdaC = invMassLcToPKPi(candLambdaC);
        deltaMass = massSigmac - massLambdaC;
        /// fill the histograms
        if (chargeSigmac == 0) {
          registry.fill(HIST("Data/hPtSigmacZero"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZero"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZero"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZero"), deltaMass, ptSigmac); // Σc0
          registry.fill(HIST("Data/hPtSoftPiSigmacZero"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZero"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZero"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtSigmacZeroPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZeroPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZeroPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZeroPlusPlus"), deltaMass, ptSigmac); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSigmacZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZero"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZero"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZero"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZero"), deltaMass, ptLambdaC); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        } else {                                                                                     /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmacZeroPlusPlus.cxx
          registry.fill(HIST("Data/hPtSigmacPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacPlusPlus"), deltaMass, ptSigmac); // Σc++
          registry.fill(HIST("Data/hPtSoftPiSigmacPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtSigmacZeroPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZeroPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZeroPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZeroPlusPlus"), deltaMass, ptSigmac); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSigmacZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        }
      } /// end candidate Λc+ → pK-π+ (and charge conjugate)
      /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
        massSigmac = invMassScRecoLcToPiKP(candSigmac);
        massLambdaC = invMassLcToPiKP(candLambdaC);
        deltaMass = massSigmac - massLambdaC;
        /// fill the histograms
        if (chargeSigmac == 0) {
          registry.fill(HIST("Data/hPtSigmacZero"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZero"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZero"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZero"), deltaMass, ptSigmac); // Σc0
          registry.fill(HIST("Data/hPtSoftPiSigmacZero"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZero"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZero"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtSigmacZeroPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZeroPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZeroPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZeroPlusPlus"), deltaMass, ptSigmac); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSigmacZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZero"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZero"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZero"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZero"), deltaMass, ptLambdaC); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        } else {                                                                                     /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmacZeroPlusPlus.cxx
          registry.fill(HIST("Data/hPtSigmacPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacPlusPlus"), deltaMass, ptSigmac); // Σc++
          registry.fill(HIST("Data/hPtSoftPiSigmacPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtSigmacZeroPlusPlus"), ptSigmac);
          registry.fill(HIST("Data/hEtaSigmacZeroPlusPlus"), etaSigmac);
          registry.fill(HIST("Data/hPhiSigmacZeroPlusPlus"), phiSigmac);
          registry.fill(HIST("Data/hDeltaMassSigmacZeroPlusPlus"), deltaMass, ptSigmac); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiSigmacZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiSigmacZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiSigmacZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLambdaCFromSigmacZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromSigmacZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromSigmacZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmacZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        }
      } /// end candidate Λc+ → π+K-p (and charge conjugate)
    }   /// end loop over the candidate Σc0,++
  };    /// end process
  PROCESS_SWITCH(HfTaskSigmac, process, "Process data", true);

  /// @brief process function to fill the histograms needed in analysis (MC)
  /// @param candidatesSigmac are the reconstructed candidate Σc0,++ with MC info
  /// @param particlesMc are the generated particles with flags wheter they are Σc0,++ or not
  /// @param
  void processMc(const soa::Join<aod::HfCandSigmac, aod::HfCandSigmacMcRec>& candidatesSigmac,
                 soa::Join<aod::McParticles, aod::HfCandSigmacMcGen> const& particlesMc,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMcLambdaC,
                 soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&, const aod::BigTracksMC&)
  {

    /// MC generated particles
    for (auto& particle : particlesMc) {

      /// reject immediately particles different from Σc0,++
      bool isSigmacZeroGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::Sigmac0ToPKPiPi));
      bool isSigmacPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::SigmacplusplusToPKPiPi));
      if (!isSigmacZeroGen && !isSigmacPlusPlusGen)
        continue;

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
      if (yCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > yCandMax) {
        continue;
      }

      /// Get the kinematic information of Σc0,++ and the daughters
      /// Get information about origin (prompt, non-prompt)
      /// Get information about decay Λc+ channel (direct, resonant)
      double ptGenSigmac(particle.pt()), etaGenSigmac(particle.eta()), phiGenSigmac(particle.phi());
      auto arrayDaughtersIds = particle.daughtersIds();
      if (arrayDaughtersIds.size() != 2) {
        /// This should never happen
        LOG(fatal) << "generated Σc0,++ has a number of daughter particles different than 2";
        continue;
      }
      double ptGenLambdaC(-1.), ptGenSoftPi(-1.);
      double etaGenLambdaC(-1.), etaGenSoftPi(-1.);
      double phiGenLambdaC(-1.), phiGenSoftPi(-1.);
      int origin = -1;
      int8_t channel = -1;
      if (std::abs(arrayDaughtersIds[0]) == pdg::Code::kLambdaCPlus) {
        /// daughter 0 is the Λc+, daughter 1 the soft π
        auto daugLambdaC = particlesMcLambdaC.rawIteratorAt(arrayDaughtersIds[0]);
        auto daugSoftPi = particlesMc.rawIteratorAt(arrayDaughtersIds[1]);
        ptGenLambdaC = daugLambdaC.pt();
        etaGenLambdaC = daugLambdaC.eta();
        phiGenLambdaC = daugLambdaC.phi();
        origin = daugLambdaC.originMcGen();
        channel = daugLambdaC.flagMcDecayChanGen();
        ptGenSoftPi = daugSoftPi.pt();
        etaGenSoftPi = daugSoftPi.eta();
        phiGenSoftPi = daugSoftPi.phi();
      } else if (std::abs(arrayDaughtersIds[0]) == kPiPlus) {
        /// daughter 0 is the soft π, daughter 1 the Λc+
        auto daugLambdaC = particlesMcLambdaC.rawIteratorAt(arrayDaughtersIds[1]);
        auto daugSoftPi = particlesMc.rawIteratorAt(arrayDaughtersIds[0]);
        ptGenLambdaC = daugLambdaC.pt();
        etaGenLambdaC = daugLambdaC.eta();
        phiGenLambdaC = daugLambdaC.phi();
        origin = daugLambdaC.originMcGen();
        channel = daugLambdaC.flagMcDecayChanGen();
        ptGenSoftPi = daugSoftPi.pt();
        etaGenSoftPi = daugSoftPi.eta();
        phiGenSoftPi = daugSoftPi.phi();
      }

      /// Fill histograms
      if (isSigmacZeroGen) {
        /// Generated Σc0 and Λc+ ← Σc0 signals
        registry.fill(HIST("MC/generated/hPtGenSigmacZeroSig"), ptGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSigmacZeroSig"), etaGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSigmacZeroSig"), phiGenSigmac, (double)origin, (double)channel); /// Generated Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSigmacZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmacZeroSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmacZeroSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmacZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmacZeroSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmacZeroSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0 signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSigmacZeroPlusPlusSig"), etaGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSigmacZeroPlusPlusSig"), phiGenSigmac, (double)origin, (double)channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmacZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmacZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmacZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmacZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0,++ signal
      } else if (isSigmacPlusPlusGen) {
        /// Generated Σc++ and Λc+ ← Σc++ signals
        registry.fill(HIST("MC/generated/hPtGenSigmacPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSigmacPlusPlusSig"), etaGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSigmacPlusPlusSig"), phiGenSigmac, (double)origin, (double)channel); /// Generated Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSigmacPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmacPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmacPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmacPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmacPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmacPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc++ signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSigmacZeroPlusPlusSig"), etaGenSigmac, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSigmacZeroPlusPlusSig"), phiGenSigmac, (double)origin, (double)channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmacZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmacZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmacZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmacZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0,++ signal
      }

    } /// end loop over generated particles

    /// reconstructed Σc0,++ matched to MC
    for (auto& candSigmac : candidatesSigmac) {

      /// Candidate selected as Σc0 and/or Σc++
      if (!(candSigmac.hfflag() & 1 << aod::hf_cand_sc::DecayType::Sigmac0ToPKPiPi) && !(candSigmac.hfflag() & 1 << aod::hf_cand_sc::DecayType::SigmacplusplusToPKPiPi)) {
        continue;
      }
      /// rapidity selection on Σc0,++
      if (yCandMax >= 0. && std::abs(ySc0(candSigmac)) > yCandMax && std::abs(yScPlusPlus(candSigmac)) > yCandMax) {
        continue;
      }

      /// electric charge
      const int chargeSigmac = candSigmac.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the Σc0
      /// and understand which mass hypotheses are possible
      const auto& candLambdaC = candSigmac.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candLambdaC, candSigmac);

      candLambdaC.flagMcDecayChanRec();

      /// Reconstructed Σc0 signal
      if (std::abs(candSigmac.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::Sigmac0ToPKPiPi && (chargeSigmac == 0)) {
        // Get the corresponding MC particle for Sigmac, found as the mother of the soft pion
        auto indexMcSigmacRec = RecoDecay::getMother(particlesMc, candSigmac.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmacMcGen>>(), pdg::Code::kSigmac0, true);
        auto particleSigmac = particlesMc.rawIteratorAt(indexMcSigmacRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = array{candLambdaC.prong0_as<aod::BigTracksMC>(), candLambdaC.prong1_as<aod::BigTracksMC>(), candLambdaC.prong2_as<aod::BigTracksMC>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSigmac.prong1_as<aod::BigTracksMC>().mcParticle();

        // const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmac.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmac.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSigmac(-1.), massLambdaC(-1.), deltaMass(-1.);
        double ptSigmac(candSigmac.pt()), ptLambdaC(candLambdaC.pt());
        double etaSigmac(candSigmac.eta()), etaLambdaC(candLambdaC.eta());
        double phiSigmac(candSigmac.phi()), phiLambdaC(candLambdaC.phi());
        double ptSoftPi(candSigmac.prong1().pt()), etaSoftPi(candSigmac.prong1().eta()), phiSoftPi(candSigmac.prong1().phi());
        double ptGenSigmac(particleSigmac.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
        int origin = candSigmac.originMcRec();
        auto channel = candLambdaC.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kProton) {
          massSigmac = invMassScRecoLcToPKPi(candSigmac);
          massLambdaC = invMassLcToPKPi(candLambdaC);
          deltaMass = massSigmac - massLambdaC;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSig"), deltaMass, ptSigmac);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSigPrompt"), deltaMass, ptSigmac);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSigNonPrompt"), deltaMass, ptSigmac); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSig"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigPrompt"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigNonPrompt"), deltaMass, ptLambdaC); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSig"), deltaMass, ptSigmac);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigPrompt"), deltaMass, ptSigmac);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmac); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSig"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kPiPlus) {
          massSigmac = invMassScRecoLcToPiKP(candSigmac);
          massLambdaC = invMassLcToPiKP(candLambdaC);
          deltaMass = massSigmac - massLambdaC;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
      }   /// end reconstructed Σc0 signal

      /// Reconstructed Σc++ signal
      else if (std::abs(candSigmac.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmacplusplusToPKPiPi && (std::abs(chargeSigmac) == 2)) {
        // Get the corresponding MC particle for Sigmac, found as the mother of the soft pion
        auto indexMcSigmacRec = RecoDecay::getMother(particlesMc, candSigmac.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmacMcGen>>(), pdg::Code::kSigmac0, true);
        auto particleSigmac = particlesMc.rawIteratorAt(indexMcSigmacRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = array{candLambdaC.prong0_as<aod::BigTracksMC>(), candLambdaC.prong1_as<aod::BigTracksMC>(), candLambdaC.prong2_as<aod::BigTracksMC>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSigmac.prong1_as<aod::BigTracksMC>().mcParticle();

        // const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmac.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmac.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSigmac(-1.), massLambdaC(-1.), deltaMass(-1.);
        double ptSigmac(candSigmac.pt()), ptLambdaC(candLambdaC.pt());
        double etaSigmac(candSigmac.eta()), etaLambdaC(candLambdaC.eta());
        double phiSigmac(candSigmac.phi()), phiLambdaC(candLambdaC.phi());
        double ptSoftPi(candSigmac.prong1().pt()), etaSoftPi(candSigmac.prong1().eta()), phiSoftPi(candSigmac.prong1().phi());
        double ptGenSigmac(particleSigmac.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
        int origin = candSigmac.originMcRec();
        auto channel = candLambdaC.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kProton) {
          massSigmac = invMassScRecoLcToPKPi(candSigmac);
          massLambdaC = invMassLcToPKPi(candLambdaC);
          deltaMass = massSigmac - massLambdaC;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kPiPlus) {
          massSigmac = invMassScRecoLcToPiKP(candSigmac);
          massLambdaC = invMassLcToPiKP(candLambdaC);
          deltaMass = massSigmac - massLambdaC;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacPlusPlusSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSigmacZeroPlusPlusSig"), ptSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSigmacZeroPlusPlusSig"), ptGenSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSigmacZeroPlusPlusSig"), etaSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSigmacZeroPlusPlusSig"), phiSigmac, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSig"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigPrompt"), deltaMass, ptSigmac, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmac, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmacZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmacZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmacZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmacZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmacZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmacZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmacZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmacZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmacZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
      }   /// end reconstructed Σc++ signal

    } /// end loop on reconstructed Σc0,++

  }; /// end processMc
  PROCESS_SWITCH(HfTaskSigmac, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSigmac>(cfgc)};
}