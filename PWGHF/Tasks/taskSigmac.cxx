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

/// \file taskSc.cxx
/// \brief Task for Σc0,++ → Λc+(→pK-π+) π-,+ analysis
/// \note Σc0,++ candidates built in O2Physics/PWGHF/TableProducer/HFCandidateCreatorScZeroPlusPlus.cxx
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

struct HfTaskSc {

  /// One value of rapidity only
  /// Remember that in Run2 the distinction among GenLimAcc, GenAccMother, GenAcc was done, where:
  ///  - GenLimAcc: Sc in |y|<0.5
  ///  - GenAccMother: Sc in the y-range of the reconstruction ("fiducial acceptance")
  ///  - GenAcc: Sc and Lc in fiducial acceptance, L daughters in acceptance
  /// Properly normalize your results to provide a cross section
  /// OR
  /// consider the new parametrization of the fiducial acceptance (to be seen for reco signal in MC)
  Configurable<float> yCandMax{"yCandMax", -1, "Sc rapidity"};

  /// analysis histograms
  HistogramRegistry registry{
    "registry",
    {/// Σc0
     {"Data/hPtScZero", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaScZero", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiScZero", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassScZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiScZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiScZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiScZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                      /// Σc++
     {"Data/hPtScPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaScPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiScPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassScPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                            /// Σc0,++
     {"Data/hPtScZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaScZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiScZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassScZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiScZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
     {"Data/hEtaSoftPiScZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                  /// soft π
     {"Data/hPhiSoftPiScZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},           /// soft π
                                                                                                                                                                                                    /// Λc+ ← Σc0
     {"Data/hPtLambdaCFromScZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromScZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromScZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromScZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc++
     {"Data/hPtLambdaCFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc0,++
     {"Data/hPtLambdaCFromScZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLambdaCFromScZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLambdaCFromScZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2 * M_PI}}}},
     {"Data/hDeltaMassLambdaCFromScZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}};

  /// @brief init function, to define the additional analysis histograms
  /// @param
  void init(InitContext&)
  {
    if (doprocessMc) {
      /////////////////////
      ///   Generated   ///
      /////////////////////
      /// Generated Σc0 signal
      registry.add("MC/generated/hPtGenScZeroSig", "#Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenScZeroSig", "#Sigma_{c}^{0} generated signal; #eta^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenScZeroSig", "#Sigma_{c}^{0} generated signal; #varphi^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Σc++ signal
      registry.add("MC/generated/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #eta^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #varphi^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Σc0,++ signal
      registry.add("MC/generated/hPtGenScZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenScZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenScZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/generated/hEtaGenSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                  /// soft π
      registry.add("MC/generated/hPhiGenSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});           /// soft π
      /// Generated Λc+ ← Σc0 signal
      registry.add("MC/generated/hPtGenLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc++ signal
      registry.add("MC/generated/hPtGenLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      /// Generated Λc+ ← Σc0,++ signal
      registry.add("MC/generated/hPtGenLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});

      /////////////////////////
      ///   Reconstructed   ///
      /////////////////////////
      /// Reconstructed Σc0 signal
      registry.add("MC/reconstructed/hPtScZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenScZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaScZeroSig", "#Sigma_{c}^{0} reconstructed signal; #eta(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiScZeroSig", "#Sigma_{c}^{0} reconstructed signal; #varphi(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroSigPrompt", "#Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroSigNonPrompt", "#Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiScZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Σc++ signal
      registry.add("MC/reconstructed/hPtScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt", "#Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt", "#Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Σc0,++ signal
      registry.add("MC/reconstructed/hPtScZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenScZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaScZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #eta(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiScZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #varphi(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroPlusPlusSigPrompt", "#Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScZeroPlusPlusSigNonPrompt", "#Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiScZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                    /// soft π
      /// Reconstructed Λc+ ← Σc0 signal
      registry.add("MC/reconstructed/hPtLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc++ signal
      registry.add("MC/reconstructed/hPtLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc0,++ signal
      registry.add("MC/reconstructed/hPtLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2 * M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
    }
  }; /// end init

  /// @brief Function to determine if the reconstructed candidate Σc0,++ decays into Λc+ → pK-π+, Λc+ → π+K-p or both
  /// @tparam L template for LambdaC daughter of Sc candidate
  /// @tparam S template for Sc candidate
  /// @param candidateLc LambdaC daughter of Sc candidate
  /// @param candSc Sc candidate
  /// @return 0: none; 1: only Λc+ → pK-π+ possible; 2: Λc+ → π+K-p possible; 3: both possible
  template <typename L, typename S>
  int isDecayToPKPiToPiKP(L& candidateLc, S& candSc)
  {
    int channel = 0;
    if ((candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG()) {
      // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      channel += 1;
    }
    if ((candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG()) {
      // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      channel += 2;
    }
    return channel; /// 0: none; 1: pK-π+ only; 2: π+K-p only; 3: both possible
  }

  using RecoLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;
  /// @brief process function to fill the histograms needed in analysis (data)
  /// @param candidatesSc are the reconstructed candidate Σc0,++
  /// @param
  void process(const aod::HfCandSc& candidatesSc,
               const RecoLc&, const aod::BigTracksExtended&)
  {

    /// loop over the candidate Σc0,++
    for (auto& candSc : candidatesSc) {

      const int chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the candidate Σc0,++
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<RecoLc>();
      // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candidateLc, candSc);
      double massSc(-1.), massLambdaC(-1.), deltaMass(-1.);
      double ptSc(candSc.pt()), ptLambdaC(candidateLc.pt());
      double etaSc(candSc.eta()), etaLambdaC(candidateLc.eta());
      double phiSc(candSc.phi()), phiLambdaC(candidateLc.phi());
      double ptSoftPi(candSc.prong1_as<aod::BigTracksExtended>().pt()), etaSoftPi(candSc.prong1_as<aod::BigTracksExtended>().eta()), phiSoftPi(candSc.prong1_as<aod::BigTracksExtended>().phi());
      /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
        massSc = invMassScRecoLcToPKPi(candSc, candidateLc);
        massLambdaC = invMassLcToPKPi(candidateLc);
        deltaMass = massSc - massLambdaC;
        /// fill the histograms
        if (chargeSc == 0) {
          registry.fill(HIST("Data/hPtScZero"), ptSc);
          registry.fill(HIST("Data/hEtaScZero"), etaSc);
          registry.fill(HIST("Data/hPhiScZero"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZero"), deltaMass, ptSc); // Σc0
          registry.fill(HIST("Data/hPtSoftPiScZero"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZero"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZero"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtScZeroPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScZeroPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScZeroPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZeroPlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiScZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromScZero"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZero"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZero"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZero"), deltaMass, ptLambdaC); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLambdaCFromScZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        } else {                                                                                     /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorScZeroPlusPlus.cxx
          registry.fill(HIST("Data/hPtScPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScPlusPlus"), deltaMass, ptSc); // Σc++
          registry.fill(HIST("Data/hPtSoftPiScPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtScZeroPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScZeroPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScZeroPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZeroPlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiScZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromScPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLambdaCFromScZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        }
      } /// end candidate Λc+ → pK-π+ (and charge conjugate)
      /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
        massSc = invMassScRecoLcToPiKP(candSc, candidateLc);
        massLambdaC = invMassLcToPiKP(candidateLc);
        deltaMass = massSc - massLambdaC;
        /// fill the histograms
        if (chargeSc == 0) {
          registry.fill(HIST("Data/hPtScZero"), ptSc);
          registry.fill(HIST("Data/hEtaScZero"), etaSc);
          registry.fill(HIST("Data/hPhiScZero"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZero"), deltaMass, ptSc); // Σc0
          registry.fill(HIST("Data/hPtSoftPiScZero"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZero"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZero"), phiSoftPi); // π ← Σc0
          registry.fill(HIST("Data/hPtScZeroPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScZeroPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScZeroPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZeroPlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiScZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromScZero"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZero"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZero"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZero"), deltaMass, ptLambdaC); // Λc+ ← Σc0
          registry.fill(HIST("Data/hPtLambdaCFromScZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        } else {                                                                                     /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorScZeroPlusPlus.cxx
          registry.fill(HIST("Data/hPtScPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScPlusPlus"), deltaMass, ptSc); // Σc++
          registry.fill(HIST("Data/hPtSoftPiScPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScPlusPlus"), phiSoftPi); // π ← Σc++
          registry.fill(HIST("Data/hPtScZeroPlusPlus"), ptSc);
          registry.fill(HIST("Data/hEtaScZeroPlusPlus"), etaSc);
          registry.fill(HIST("Data/hPhiScZeroPlusPlus"), phiSc);
          registry.fill(HIST("Data/hDeltaMassScZeroPlusPlus"), deltaMass, ptSc); // Σc0,++
          registry.fill(HIST("Data/hPtSoftPiScZeroPlusPlus"), ptSoftPi);
          registry.fill(HIST("Data/hEtaSoftPiScZeroPlusPlus"), etaSoftPi);
          registry.fill(HIST("Data/hPhiSoftPiScZeroPlusPlus"), phiSoftPi); // π ← Σc0,++
          registry.fill(HIST("Data/hPtLambdaCFromScPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc++
          registry.fill(HIST("Data/hPtLambdaCFromScZeroPlusPlus"), ptLambdaC);
          registry.fill(HIST("Data/hEtaLambdaCFromScZeroPlusPlus"), etaLambdaC);
          registry.fill(HIST("Data/hPhiLambdaCFromScZeroPlusPlus"), phiLambdaC);
          registry.fill(HIST("Data/hDeltaMassLambdaCFromScZeroPlusPlus"), deltaMass, ptLambdaC); // Λc+ ← Σc0,++
        }
      } /// end candidate Λc+ → π+K-p (and charge conjugate)
    }   /// end loop over the candidate Σc0,++
  };    /// end process

  /// @brief process function to fill the histograms needed in analysis (MC)
  /// @param candidatesSc are the reconstructed candidate Σc0,++ with MC info
  /// @param particlesMc are the generated particles with flags wheter they are Σc0,++ or not
  /// @param
  void processMc(const soa::Join<aod::HfCandSc, aod::HfCandScMcRec>& candidatesSc,
                 aod::McParticles const& particlesMc,
                 soa::Join<aod::McParticles, aod::HfCandScMcGen> const& particlesMcSc,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMcLc,
                 soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&, const aod::BigTracksMC&)
  {

    /// MC generated particles
    for (auto& particle : particlesMcSc) {

      /// reject immediately particles different from Σc0,++
      bool isScZeroGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::Sc0ToPKPiPi));
      bool isScPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::ScplusplusToPKPiPi));
      if (!isScZeroGen && !isScPlusPlusGen)
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
      double ptGenSc(particle.pt()), etaGenSc(particle.eta()), phiGenSc(particle.phi());
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
        auto daugLambdaC = particlesMcLc.rawIteratorAt(arrayDaughtersIds[0]);
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
        auto daugLambdaC = particlesMcLc.rawIteratorAt(arrayDaughtersIds[1]);
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
      if (isScZeroGen) {
        /// Generated Σc0 and Λc+ ← Σc0 signals
        registry.fill(HIST("MC/generated/hPtGenScZeroSig"), ptGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenScZeroSig"), etaGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenScZeroSig"), phiGenSc, (double)origin, (double)channel); /// Generated Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiScZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiScZeroSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiScZeroSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0 signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromScZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromScZeroSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromScZeroSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0 signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenScZeroPlusPlusSig"), etaGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenScZeroPlusPlusSig"), phiGenSc, (double)origin, (double)channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiScZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiScZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromScZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromScZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0,++ signal
      } else if (isScPlusPlusGen) {
        /// Generated Σc++ and Λc+ ← Σc++ signals
        registry.fill(HIST("MC/generated/hPtGenScPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenScPlusPlusSig"), etaGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenScPlusPlusSig"), phiGenSc, (double)origin, (double)channel); /// Generated Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiScPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiScPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromScPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromScPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromScPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc++ signal
        /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
        registry.fill(HIST("MC/generated/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenScZeroPlusPlusSig"), etaGenSc, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenScZeroPlusPlusSig"), phiGenSc, (double)origin, (double)channel); /// Generated Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenSoftPiScZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenSoftPiScZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double)channel); /// Generated π ← Σc0,++ signal
        registry.fill(HIST("MC/generated/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hEtaGenLambdaCFromScZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double)channel);
        registry.fill(HIST("MC/generated/hPhiGenLambdaCFromScZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double)channel); /// Generated Λc+ ← Σc0,++ signal
      }

    } /// end loop over generated particles

    /// reconstructed Σc0,++ matched to MC
    for (auto& candSc : candidatesSc) {

      /// Candidate selected as Σc0 and/or Σc++
      if (!(candSc.hfflag() & 1 << aod::hf_cand_sc::DecayType::Sc0ToPKPiPi) && !(candSc.hfflag() & 1 << aod::hf_cand_sc::DecayType::ScplusplusToPKPiPi)) {
        continue;
      }
      /// rapidity selection on Σc0,++
      if (yCandMax >= 0. && std::abs(ySc0(candSc)) > yCandMax && std::abs(yScPlusPlus(candSc)) > yCandMax) {
        continue;
      }

      /// electric charge
      const int chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the Σc0
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candidateLc, candSc);

      //candidateLc.flagMcDecayChanRec();

      /// Reconstructed Σc0 signal
      if (std::abs(candSc.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::Sc0ToPKPiPi && (chargeSc == 0)) {
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        auto indexMcScRec = RecoDecay::getMother(particlesMc, candSc.prong1_as<aod::BigTracksMC>().mcParticle(), pdg::Code::kSigmac0, true);
        auto particleSc = particlesMc.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = array{candidateLc.prong0_as<aod::BigTracksMC>(), candidateLc.prong1_as<aod::BigTracksMC>(), candidateLc.prong2_as<aod::BigTracksMC>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSc.prong1_as<aod::BigTracksMC>().mcParticle();

        // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSc(-1.), massLambdaC(-1.), deltaMass(-1.);
        double ptSc(candSc.pt()), ptLambdaC(candidateLc.pt());
        double etaSc(candSc.eta()), etaLambdaC(candidateLc.eta());
        double phiSc(candSc.phi()), phiLambdaC(candidateLc.phi());
        double ptSoftPi(candSc.prong1_as<aod::BigTracksMC>().pt()), etaSoftPi(candSc.prong1_as<aod::BigTracksMC>().eta()), phiSoftPi(candSc.prong1_as<aod::BigTracksMC>().phi());
        double ptGenSc(particleSc.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
        int origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kProton) {
          massSc = invMassScRecoLcToPKPi(candSc, candidateLc);
          massLambdaC = invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLambdaC;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSig"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSigPrompt"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSigNonPrompt"), deltaMass, ptSc); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSig"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigPrompt"), deltaMass, ptLambdaC);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigNonPrompt"), deltaMass, ptLambdaC); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSig"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigPrompt"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigNonPrompt"), deltaMass, ptSc); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = invMassScRecoLcToPiKP(candSc, candidateLc);
          massLambdaC = invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLambdaC;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
      }   /// end reconstructed Σc0 signal

      /// Reconstructed Σc++ signal
      else if (std::abs(candSc.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::ScplusplusToPKPiPi && (std::abs(chargeSc) == 2)) {
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        auto indexMcScRec = RecoDecay::getMother(particlesMc, candSc.prong1_as<aod::BigTracksMC>().mcParticle(), pdg::Code::kSigmacPlusPlus, true);
        auto particleSc = particlesMc.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = array{candidateLc.prong0_as<aod::BigTracksMC>(), candidateLc.prong1_as<aod::BigTracksMC>(), candidateLc.prong2_as<aod::BigTracksMC>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
        // Get the corresponding MC particle for soft pion
        auto particleSoftPi = candSc.prong1_as<aod::BigTracksMC>().mcParticle();

        // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
        // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
        double massSc(-1.), massLambdaC(-1.), deltaMass(-1.);
        double ptSc(candSc.pt()), ptLambdaC(candidateLc.pt());
        double etaSc(candSc.eta()), etaLambdaC(candidateLc.eta());
        double phiSc(candSc.phi()), phiLambdaC(candidateLc.phi());
        double ptSoftPi(candSc.prong1_as<aod::BigTracksMC>().pt()), etaSoftPi(candSc.prong1_as<aod::BigTracksMC>().eta()), phiSoftPi(candSc.prong1_as<aod::BigTracksMC>().phi());
        double ptGenSc(particleSc.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
        int origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kProton) {
          massSc = invMassScRecoLcToPKPi(candSc, candidateLc);
          massLambdaC = invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLambdaC;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = invMassScRecoLcToPiKP(candSc, candidateLc);
          massLambdaC = invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLambdaC;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtScZeroPlusPlusSig"), ptSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScZeroPlusPlusSig"), ptGenSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaScZeroPlusPlusSig"), etaSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiScZeroPlusPlusSig"), phiSc, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSig"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigPrompt"), deltaMass, ptSc, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScZeroPlusPlusSigNonPrompt"), deltaMass, ptSc, (double)channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScZeroPlusPlusSig"), ptSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScZeroPlusPlusSig"), etaSoftPi, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScZeroPlusPlusSig"), phiSoftPi, (double)origin, (double)channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLambdaCFromScZeroPlusPlusSig"), ptLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromScZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromScZeroPlusPlusSig"), etaLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromScZeroPlusPlusSig"), phiLambdaC, (double)origin, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSig"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double)channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromScZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double)channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
      }   /// end reconstructed Σc++ signal

    } /// end loop on reconstructed Σc0,++

  }; /// end processMc
  PROCESS_SWITCH(HfTaskSc, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSc>(cfgc)};
}