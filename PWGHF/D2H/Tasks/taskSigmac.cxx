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

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

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
  Configurable<float> yCandMax{"yCandMax", -1, "Sc rapidity"};

  HfHelper hfHelper;

  /// analysis histograms
  HistogramRegistry registry{
    "registry",
    {/// Σc0
     {"Data/hPtSc0", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSc0", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSc0", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassSc0", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiSc0", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                             /// Σc++
     {"Data/hPtScPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaScPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiScPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassScPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiScPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                                      /// Σc0,++
     {"Data/hPtSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassSc0PlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     {"Data/hPtSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}},     /// soft π
     {"Data/hEtaSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},                      /// soft π
     {"Data/hPhiSoftPiSc0PlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}}, /// soft π
                                                                                                                                                                                                           /// Λc+ ← Σc0
     {"Data/hPtLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassLcFromSc0", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc++
     {"Data/hPtLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassLcFromScPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
     /// Λc+ ← Σc0,++
     {"Data/hPtLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
     {"Data/hEtaLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
     {"Data/hPhiLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, constants::math::TwoPI}}}},
     {"Data/hDeltaMassLcFromSc0PlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}};

  using RecoLc = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;

  /// @brief init function, to define the additional analysis histograms
  /// @param
  void init(InitContext&)
  {
    if (doprocessMc) {
      /////////////////////
      ///   Generated   ///
      /////////////////////
      /// Generated Σc0 signal
      registry.add("MC/generated/hPtGenSc0Sig", "#Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSc0Sig", "#Sigma_{c}^{0} generated signal; #eta^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSc0Sig", "#Sigma_{c}^{0} generated signal; #varphi^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                       /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});  /// soft π
      /// Generated Σc++ signal
      registry.add("MC/generated/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #eta^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenScPlusPlusSig", "#Sigma_{c}^{++} generated signal; #varphi^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      registry.add("MC/generated/hEtaGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                       /// soft π
      registry.add("MC/generated/hPhiGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});  /// soft π
      /// Generated Σc0,++ signal
      registry.add("MC/generated/hPtGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hEtaGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPhiGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/generated/hPtGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      registry.add("MC/generated/hEtaGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                       /// soft π
      registry.add("MC/generated/hPhiGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});  /// soft π
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
      registry.add("MC/reconstructed/hDeltaMassSc0Sig", "#Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0SigPrompt", "#Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0SigNonPrompt", "#Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSc0Sig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Σc++ signal
      registry.add("MC/reconstructed/hPtScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt", "#Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt", "#Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiScPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Σc0,++ signal
      registry.add("MC/reconstructed/hPtSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #eta(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #varphi(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt", "#Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt", "#Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});          /// soft π
      registry.add("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
      registry.add("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});                           /// soft π
      registry.add("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});      /// soft π
      /// Reconstructed Λc+ ← Σc0 signal
      registry.add("MC/reconstructed/hPtLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0Sig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc++ signal
      registry.add("MC/reconstructed/hPtLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      /// Reconstructed Λc+ ← Σc0,++ signal
      registry.add("MC/reconstructed/hPtLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hEtaLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hPhiLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, constants::math::TwoPI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
      registry.add("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
    }
  }; /// end init

  /// @brief Function to determine if the reconstructed candidate Σc0,++ decays into Λc+ → pK-π+, Λc+ → π+K-p or both
  /// @tparam L template for Lc daughter of Sc candidate
  /// @tparam S template for Sc candidate
  /// @param candidateLc Lc daughter of Sc candidate
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

  /// @brief process function to fill the histograms needed in analysis (data)
  /// @param candidatesSc are the reconstructed candidate Σc0,++
  /// @param
  void process(aod::HfCandSc const& candidatesSc,
               RecoLc const&,
               aod::Tracks const&)
  {

    /// loop over the candidate Σc0,++
    for (const auto& candSc : candidatesSc) {

      const int8_t chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the candidate Σc0,++
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<RecoLc>();
      // const int iscandidateLcpKpi = (candidateLc.isSelLcToPKPi() >= 1) && candSc.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
      // const int iscandidateLcpiKp = (candidateLc.isSelLcToPiKP() >= 1) && candSc.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candidateLc, candSc);
      double massSc(-1.), massLc(-1.), deltaMass(-1.);
      double ptSc(candSc.pt()), ptLc(candidateLc.pt());
      double etaSc(candSc.eta()), etaLc(candidateLc.eta());
      double phiSc(candSc.phi()), phiLc(candidateLc.phi());
      double ptSoftPi(candSc.prong1().pt()), etaSoftPi(candSc.prong1().eta()), phiSoftPi(candSc.prong1().phi());
      /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
        massSc = hfHelper.invMassScRecoLcToPKPi(candSc, candidateLc);
        massLc = hfHelper.invMassLcToPKPi(candidateLc);
        deltaMass = massSc - massLc;
        /// fill the histograms
        if (chargeSc == 0) {
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
      } /// end candidate Λc+ → pK-π+ (and charge conjugate)
      /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
      if (isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
        massSc = hfHelper.invMassScRecoLcToPiKP(candSc, candidateLc);
        massLc = hfHelper.invMassLcToPiKP(candidateLc);
        deltaMass = massSc - massLc;
        /// fill the histograms
        if (chargeSc == 0) {
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
      } /// end candidate Λc+ → π+K-p (and charge conjugate)
    }   /// end loop over the candidate Σc0,++
  };    /// end process

  /// @brief process function to fill the histograms needed in analysis (MC)
  /// @param candidatesSc are the reconstructed candidate Σc0,++ with MC info
  /// @param mcParticles are the generated particles with flags wheter they are Σc0,++ or not
  /// @param
  void processMc(soa::Join<aod::HfCandSc, aod::HfCandScMcRec> const& candidatesSc,
                 aod::McParticles const& mcParticles,
                 soa::Join<aod::McParticles, aod::HfCandScMcGen> const& mcParticlesSc,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticlesLc,
                 soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&,
                 aod::TracksWMc const&)
  {

    /// MC generated particles
    for (const auto& particle : mcParticlesSc) {

      /// reject immediately particles different from Σc0,++
      bool isSc0Gen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi));
      bool isScPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi));
      if (!isSc0Gen && !isScPlusPlusGen)
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
      if (yCandMax >= 0. && std::abs(RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassSigmaC0)) > yCandMax) {
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
      double ptGenLc(-1.), ptGenSoftPi(-1.);
      double etaGenLc(-1.), etaGenSoftPi(-1.);
      double phiGenLc(-1.), phiGenSoftPi(-1.);
      int origin = -1;
      int8_t channel = -1;
      if (std::abs(arrayDaughtersIds[0]) == pdg::Code::kLambdaCPlus) {
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
      } else if (std::abs(arrayDaughtersIds[0]) == kPiPlus) {
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

      /// Fill histograms
      if (isSc0Gen) {
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
      } else if (isScPlusPlusGen) {
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
      }

    } /// end loop over generated particles

    /// reconstructed Σc0,++ matched to MC
    for (const auto& candSc : candidatesSc) {

      /// Candidate selected as Σc0 and/or Σc++
      if (!(candSc.hfflag() & 1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi) && !(candSc.hfflag() & 1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi)) {
        continue;
      }
      /// rapidity selection on Σc0,++
      if (yCandMax >= 0. && std::abs(hfHelper.ySc0(candSc)) > yCandMax && std::abs(hfHelper.yScPlusPlus(candSc)) > yCandMax) {
        continue;
      }

      /// electric charge
      const int8_t chargeSc = candSc.charge(); // either Σc0 or Σc++

      /// get the candidate Λc+ used to build the Σc0
      /// and understand which mass hypotheses are possible
      const auto& candidateLc = candSc.prongLc_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
      const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candidateLc, candSc);

      // candidateLc.flagMcDecayChanRec();

      /// Reconstructed Σc0 signal
      if (std::abs(candSc.flagMcMatchRec()) == 1 << aod::hf_cand_sigmac::DecayType::Sc0ToPKPiPi && (chargeSc == 0)) {
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        auto indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), pdg::Code::kSigmaC0, true);
        auto particleSc = mcParticles.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = std::array{candidateLc.prong0_as<aod::TracksWMc>(), candidateLc.prong1_as<aod::TracksWMc>(), candidateLc.prong2_as<aod::TracksWMc>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersLc, pdg::Code::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
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
        int origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kProton) {
          massSc = hfHelper.invMassScRecoLcToPKPi(candSc, candidateLc);
          massLc = hfHelper.invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLc;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSc0Sig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0Sig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0Sig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0Sig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0Sig"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigPrompt"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigNonPrompt"), deltaMass, ptSc); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0Sig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0Sig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0Sig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0Sig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0Sig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0Sig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0Sig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0Sig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0Sig"), deltaMass, ptLc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt"), deltaMass, ptLc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt"), deltaMass, ptLc); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = hfHelper.invMassScRecoLcToPiKP(candSc, candidateLc);
          massLc = hfHelper.invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLc;

          /// Fill the histograms for reconstructed Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSc0Sig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0Sig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0Sig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0Sig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0Sig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0SigNonPrompt"), deltaMass, ptSc, channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0Sig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0Sig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0Sig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0Sig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0Sig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0Sig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0Sig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0Sig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0Sig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0SigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → π+K-p (and charge conjugate)
        /// end reconstructed Σc0 signal
      } else if (std::abs(candSc.flagMcMatchRec()) == 1 << aod::hf_cand_sigmac::DecayType::ScplusplusToPKPiPi && (std::abs(chargeSc) == 2)) {
        /// Reconstructed Σc++ signal
        // Get the corresponding MC particle for Sc, found as the mother of the soft pion
        auto indexMcScRec = RecoDecay::getMother(mcParticles, candSc.prong1_as<aod::TracksWMc>().mcParticle(), pdg::Code::kSigmaCPlusPlus, true);
        auto particleSc = mcParticles.rawIteratorAt(indexMcScRec);
        // Get the corresponding MC particle for Lc
        auto arrayDaughtersLc = std::array{candidateLc.prong0_as<aod::TracksWMc>(), candidateLc.prong1_as<aod::TracksWMc>(), candidateLc.prong2_as<aod::TracksWMc>()};
        int8_t sign = 0;
        int indexMcLcRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersLc, pdg::Code::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
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
        int origin = candSc.originMcRec();
        auto channel = candidateLc.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

        /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kProton) {
          massSc = hfHelper.invMassScRecoLcToPKPi(candSc, candidateLc);
          massLc = hfHelper.invMassLcToPKPi(candidateLc);
          deltaMass = massSc - massLc;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromScPlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromScPlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromScPlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromScPlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal

        } /// end candidate Λc+ → pK-π+ (and charge conjugate)
        /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
        if ((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candidateLc.prong0_as<aod::TracksWMc>().mcParticle().pdgCode()) == kPiPlus) {
          massSc = hfHelper.invMassScRecoLcToPiKP(candSc, candidateLc);
          massLc = hfHelper.invMassLcToPiKP(candidateLc);
          deltaMass = massSc - massLc;

          /// Fill the histograms for reconstructed Σc++ signal
          registry.fill(HIST("MC/reconstructed/hPtScPlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenScPlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaScPlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiScPlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassScPlusPlusSigNonPrompt"), deltaMass, ptSc, channel); // Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiScPlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiScPlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiScPlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiScPlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0 signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromScPlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromScPlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromScPlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromScPlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromScPlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0 signal

          /// Fill the histograms for reconstructed Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSc0PlusPlusSig"), ptSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSc0PlusPlusSig"), ptGenSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSc0PlusPlusSig"), etaSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSc0PlusPlusSig"), phiSc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSig"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigPrompt"), deltaMass, ptSc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassSc0PlusPlusSigNonPrompt"), deltaMass, ptSc, channel); // Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtSoftPiSc0PlusPlusSig"), ptSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSc0PlusPlusSig"), ptGenSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaSoftPiSc0PlusPlusSig"), etaSoftPi, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiSoftPiSc0PlusPlusSig"), phiSoftPi, origin, channel); // π ← Σc0,++ signal
          registry.fill(HIST("MC/reconstructed/hPtLcFromSc0PlusPlusSig"), ptLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPtGenLcFromSc0PlusPlusSig"), ptGenLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hEtaLcFromSc0PlusPlusSig"), etaLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hPhiLcFromSc0PlusPlusSig"), phiLc, origin, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSig"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigPrompt"), deltaMass, ptLc, channel);
          registry.fill(HIST("MC/reconstructed/hDeltaMassLcFromSc0PlusPlusSigNonPrompt"), deltaMass, ptLc, channel); // Λc+ ← Σc0,++ signal

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
