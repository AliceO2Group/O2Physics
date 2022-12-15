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

/// \file taskSigmaC.cxx
/// \brief Task for Σc0,++ → Λc+(→pK-π+) π-,+ analysis
/// \note Σc0,++ candidates built in O2Physics/PWGHF/TableProducer/HFCandidateCreatorSigmaCZeroPlusPlus.cxx
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

struct HfTaskSigmaC{

    /// One value of rapidity only
    /// Remember that in Run2 the distinction among GenLimAcc, GenAccMother, GenAcc was done, where:
    ///  - GenLimAcc: Sc in |y|<0.5
    ///  - GenAccMother: Sc in the y-range of the reconstruction ("fiducial acceptance")
    ///  - GenAcc: Sc and Lc in fiducial acceptance, L daughters in acceptance
    /// Properly normalize your results to provide a cross section
    /// OR
    /// consider the new parametrization of the fiducial acceptance (to be seen for reco signal in MC)
    Configurable<float> yCandMax{"yCandMax", -1, "SigmaC rapidity"};

    /// analysis histograms
    HistogramRegistry registry{
        "registry",
        {/// Σc0
        {"Data/hPtSigmaCZero", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaSigmaCZero", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiSigmaCZero", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Σc++
        {"Data/hPtSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Σc0,++
        {"Data/hPtSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{60, 0., 6.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Λc+ ← Σc0
        {"Data/hPtLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc++
        {"Data/hPtLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0,++
        {"Data/hPtLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{40, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}
    };

    /// @brief init function, to define the additional analysis histograms
    /// @param  
    void init(InitContext&) {
        if(doprocessMc) {
            /////////////////////
            ///   Generated   ///
            /////////////////////
            /// Generated Σc0 signal
            registry.add("MC/generated/hPtGenSigmaCZeroSig", "#Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenSigmaCZeroSig", "#Sigma_{c}^{0} generated signal; #eta^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenSigmaCZeroSig", "#Sigma_{c}^{0} generated signal; #varphi^{gen}(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPtGenSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hEtaGenSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hPhiGenSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Generated Σc++ signal
            registry.add("MC/generated/hPtGenSigmaCPlusPlusSig", "#Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenSigmaCPlusPlusSig", "#Sigma_{c}^{++} generated signal; #eta^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenSigmaCPlusPlusSig", "#Sigma_{c}^{++} generated signal; #varphi^{gen}(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPtGenSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hEtaGenSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hPhiGenSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Generated Σc0,++ signal
            registry.add("MC/generated/hPtGenSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPtGenSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hEtaGenSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/generated/hPhiGenSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Generated Λc+ ← Σc0 signal
            registry.add("MC/generated/hPtGenLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            /// Generated Λc+ ← Σc++ signal
            registry.add("MC/generated/hPtGenLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            /// Generated Λc+ ← Σc0,++ signal
            registry.add("MC/generated/hPtGenLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hEtaGenLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #eta^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/generated/hPhiGenLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} generated signal; #varphi^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});

            /////////////////////////
            ///   Reconstructed   ///
            /////////////////////////
            /// Reconstructed Σc0 signal
            registry.add("MC/reconstructed/hPtSigmaCZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenSigmaCZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaSigmaCZeroSig", "#Sigma_{c}^{0} reconstructed signal; #eta(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiSigmaCZeroSig", "#Sigma_{c}^{0} reconstructed signal; #varphi(#Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSig", "#Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt", "#Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt", "#Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPtGenSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hEtaSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPhiSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Reconstructed Σc++ signal
            registry.add("MC/reconstructed/hPtSigmaCPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenSigmaCPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaSigmaCPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiSigmaCPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig", "#Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt", "#Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt", "#Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPtGenSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Reconstructed Σc0,++ signal
            registry.add("MC/reconstructed/hPtSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #eta(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #varphi(#Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSig", "#Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigPrompt", "#Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigNonPrompt", "#Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPtGenSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{60, 0., 6.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hEtaSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            registry.add("MC/reconstructed/hPhiSoftPiSigmaCZeroPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}}); /// soft π
            /// Reconstructed Λc+ ← Σc0 signal
            registry.add("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            /// Reconstructed Λc+ ← Σc++ signal
            registry.add("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}, {4, -0.5, 3.5}}});
            /// Reconstructed Λc+ ← Σc0,++ signal
            registry.add("MC/reconstructed/hPtLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{p}_{T}^{gen}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{360, 0., 36.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hEtaLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{40, -2., 2.}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hPhiLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); origin (1 : prompt, 2: non-prompt); channel (0: direct  1,2,3: resonant);", {HistType::kTH3D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}, {4, -0.5, 3.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} reconstructed non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH3D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
        }
    }; /// end init

    /// @brief Function to determine if the reconstructed candidate Σc0,++ decays into Λc+ → pK-π+, Λc+ → π+K-p or both
    /// @tparam L template for LambdaC daughter of SigmaC candidate
    /// @tparam S template for SigmaC candidate
    /// @param candLambdaC LambdaC daughter of SigmaC candidate
    /// @param candSigmaC SigmaC candidate
    /// @return 0: none; 1: only Λc+ → pK-π+ possible; 2: Λc+ → π+K-p possible; 3: both possible
    template <typename L, typename S>
    int isDecayToPKPiToPiKP(L& candLambdaC, S& candSigmaC){
        int channel = 0;
        if( (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG() ){
            // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
            channel += 1;
        }
        if( (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG() ){
            // Λc+ → π+K-p and within the requested mass to build the Σc0,++
            channel += 2;
        }
        return channel; /// 0: none; 1: pK-π+ only; 2: π+K-p only; 3: both possible
    }


    /// @brief process function to fill the histograms needed in analysis (data)
    /// @param candidatesSigmaC are the reconstructed candidate Σc0,++
    /// @param 
    void process(const aod::HfCandSigmaC& candidatesSigmaC,
    soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&, const soa::Join<aod::Tracks, aod::TracksDCA>&) {

        /// loop over the candidate Σc0,++
        for(auto& candSigmaC : candidatesSigmaC) {

            const int chargeSigmaC = candSigmaC.charge();   // either Σc0 or Σc++

            /// get the candidate Λc+ used to build the candidate Σc0,++
            /// and understand which mass hypotheses are possible
            const auto& candLambdaC = candSigmaC.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
            //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
            //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
            const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candLambdaC, candSigmaC);
            double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.);
            double ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
            double etaSigmaC(candSigmaC.eta()), etaLambdaC(candLambdaC.eta());
            double phiSigmaC(candSigmaC.phi()), phiLambdaC(candLambdaC.phi());
            double ptSoftPi(candSigmaC.prong1().pt()), etaSoftPi(candSigmaC.prong1().eta()), phiSoftPi(candSigmaC.prong1().phi());
            /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
            if(isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
                massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                massLambdaC = invMassLcToPKPi(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("Data/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZero"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZero"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZero"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZero"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZero"), phiSoftPi);   // π ← Σc0
                    registry.fill(HIST("Data/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZeroPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZeroPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZeroPlusPlus"), phiSoftPi);   // π ← Σc0,++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZero"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZero"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("Data/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCPlusPlus"), phiSoftPi);   // π ← Σc++
                    registry.fill(HIST("Data/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZeroPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZeroPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZeroPlusPlus"), phiSoftPi);   // π ← Σc0,++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → pK-π+ (and charge conjugate)
            /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
            if(isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
                massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                massLambdaC = invMassLcToPiKP(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("Data/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZero"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZero"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZero"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZero"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZero"), phiSoftPi);   // π ← Σc0
                    registry.fill(HIST("Data/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZeroPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZeroPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZeroPlusPlus"), phiSoftPi);   // π ← Σc0,++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZero"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZero"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("Data/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCPlusPlus"), phiSoftPi);   // π ← Σc++
                    registry.fill(HIST("Data/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("Data/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("Data/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("Data/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("Data/hPtSoftPiSigmaCZeroPlusPlus"), ptSoftPi);
                    registry.fill(HIST("Data/hEtaSoftPiSigmaCZeroPlusPlus"), etaSoftPi);
                    registry.fill(HIST("Data/hPhiSoftPiSigmaCZeroPlusPlus"), phiSoftPi);   // π ← Σc0,++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("Data/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("Data/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("Data/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → π+K-p (and charge conjugate)
        } /// end loop over the candidate Σc0,++
    };  /// end process
    PROCESS_SWITCH(HfTaskSigmaC, process, "Process data", true);

    /// @brief process function to fill the histograms needed in analysis (MC)
    /// @param candidatesSigmaC are the reconstructed candidate Σc0,++ with MC info
    /// @param particlesMc are the generated particles with flags wheter they are Σc0,++ or not
    /// @param 
    void processMc(const soa::Join<aod::HfCandSigmaC, aod::HfCandSigmaCMcRec>& candidatesSigmaC,
    soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen> const& particlesMc,
    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMcLambdaC,
    soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const&, const aod::BigTracksMC&) {

        /// MC generated particles
        for(auto& particle : particlesMc) {

            /// reject immediately particles different from Σc0,++
            bool isSigmaCZeroGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi));
            bool isSigmaCPlusPlusGen = (std::abs(particle.flagMcMatchGen()) == (1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi));
            if(!isSigmaCZeroGen && !isSigmaCPlusPlusGen)
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
            double ptGenSigmaC(particle.pt()), etaGenSigmaC(particle.eta()), phiGenSigmaC(particle.phi());
            auto arrayDaughtersIds = particle.daughtersIds();
            if (arrayDaughtersIds.size() !=2) {
                /// This should never happen
                LOG(fatal) << "generated Σc0,++ has a number of daughter particles different than 2";
                continue;
            }
            double ptGenLambdaC(-1.), ptGenSoftPi(-1.);
            double etaGenLambdaC(-1.), etaGenSoftPi(-1.);
            double phiGenLambdaC(-1.), phiGenSoftPi(-1.);
            int origin = -1;
            int8_t channel = -1;
            if(std::abs(arrayDaughtersIds[0]) == pdg::Code::kLambdaCPlus) {
                /// daughter 0 is the Λc+, daughter 1 the soft π
                auto daugLambdaC = particlesMcLambdaC.rawIteratorAt( arrayDaughtersIds[0] );
                auto daugSoftPi = particlesMc.rawIteratorAt( arrayDaughtersIds[1] );
                ptGenLambdaC = daugLambdaC.pt();
                etaGenLambdaC = daugLambdaC.eta();
                phiGenLambdaC = daugLambdaC.phi();
                origin = daugLambdaC.originMcGen();
                channel = daugLambdaC.flagMcDecayChanGen();
                ptGenSoftPi = daugSoftPi.pt();
                etaGenSoftPi = daugSoftPi.eta();
                phiGenSoftPi = daugSoftPi.phi();
            }
            else if(std::abs(arrayDaughtersIds[0]) == kPiPlus) {
                /// daughter 0 is the soft π, daughter 1 the Λc+ 
                auto daugLambdaC = particlesMcLambdaC.rawIteratorAt( arrayDaughtersIds[1] );
                auto daugSoftPi = particlesMc.rawIteratorAt( arrayDaughtersIds[0] );
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
            if (isSigmaCZeroGen) {
                /// Generated Σc0 and Λc+ ← Σc0 signals
                registry.fill(HIST("MC/generated/hPtGenSigmaCZeroSig"), ptGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSigmaCZeroSig"), etaGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSigmaCZeroSig"), phiGenSigmaC, (double)origin, (double) channel); /// Generated Σc0 signal
                registry.fill(HIST("MC/generated/hPtGenSoftPiSigmaCZeroSig"), ptGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmaCZeroSig"), etaGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmaCZeroSig"), phiGenSoftPi, (double)origin, (double) channel); /// Generated π ← Σc0 signal
                registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmaCZeroSig"), ptGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmaCZeroSig"), etaGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmaCZeroSig"), phiGenLambdaC, (double)origin, (double) channel); /// Generated Λc+ ← Σc0 signal
                /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
                registry.fill(HIST("MC/generated/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSigmaCZeroPlusPlusSig"), etaGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSigmaCZeroPlusPlusSig"), phiGenSigmaC, (double)origin, (double) channel); /// Generated Σc0,++ signal
                registry.fill(HIST("MC/generated/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmaCZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmaCZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double) channel); /// Generated π ← Σc0,++ signal
                registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmaCZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmaCZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double) channel); /// Generated Λc+ ← Σc0,++ signal
            }
            else if (isSigmaCPlusPlusGen) {
                /// Generated Σc++ and Λc+ ← Σc++ signals
                registry.fill(HIST("MC/generated/hPtGenSigmaCPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSigmaCPlusPlusSig"), etaGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSigmaCPlusPlusSig"), phiGenSigmaC, (double)origin, (double) channel); /// Generated Σc++ signal
                registry.fill(HIST("MC/generated/hPtGenSoftPiSigmaCPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmaCPlusPlusSig"), etaGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmaCPlusPlusSig"), phiGenSoftPi, (double)origin, (double) channel); /// Generated π ← Σc++ signal
                registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmaCPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmaCPlusPlusSig"), etaGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmaCPlusPlusSig"), phiGenLambdaC, (double)origin, (double) channel); /// Generated Λc+ ← Σc++ signal
                /// Generated Σc0,++ and Λc+ ← Σc0,++ signals
                registry.fill(HIST("MC/generated/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSigmaCZeroPlusPlusSig"), etaGenSigmaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSigmaCZeroPlusPlusSig"), phiGenSigmaC, (double)origin, (double) channel); /// Generated Σc0,++ signal
                registry.fill(HIST("MC/generated/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenSoftPiSigmaCZeroPlusPlusSig"), etaGenSoftPi, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenSoftPiSigmaCZeroPlusPlusSig"), phiGenSoftPi, (double)origin, (double) channel); /// Generated π ← Σc0,++ signal
                registry.fill(HIST("MC/generated/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hEtaGenLambdaCFromSigmaCZeroPlusPlusSig"), etaGenLambdaC, (double)origin, (double) channel);
                registry.fill(HIST("MC/generated/hPhiGenLambdaCFromSigmaCZeroPlusPlusSig"), phiGenLambdaC, (double)origin, (double) channel); /// Generated Λc+ ← Σc0,++ signal
            }
                        

        } /// end loop over generated particles

        /// reconstructed Σc0,++ matched to MC
        for(auto& candSigmaC : candidatesSigmaC) {

            /// Candidate selected as Σc0 and/or Σc++
            if (!(candSigmaC.hfflag() & 1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi) && !(candSigmaC.hfflag() & 1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi)) {
              continue;
            }
            /// rapidity selection on Σc0,++
            if (yCandMax >= 0. && std::abs(ySc0(candSigmaC)) > yCandMax && std::abs(yScPlusPlus(candSigmaC)) > yCandMax) {
              continue;
            }

            /// electric charge
            const int chargeSigmaC = candSigmaC.charge();   // either Σc0 or Σc++

            /// get the candidate Λc+ used to build the Σc0
            /// and understand which mass hypotheses are possible
            const auto& candLambdaC = candSigmaC.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
            const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candLambdaC, candSigmaC);

            candLambdaC.flagMcDecayChanRec();

            /// Reconstructed Σc0 signal
            if(std::abs(candSigmaC.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi && (chargeSigmaC == 0)) {
                // Get the corresponding MC particle for SigmaC, found as the mother of the soft pion
                auto indexMcSigmaCRec = RecoDecay::getMother(particlesMc, candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen>>(), pdg::Code::kSigmaC0, true);
                auto particleSigmaC = particlesMc.rawIteratorAt(indexMcSigmaCRec);
                // Get the corresponding MC particle for Lc
                auto arrayDaughtersLc = array{candLambdaC.prong0_as<aod::BigTracksMC>(), candLambdaC.prong1_as<aod::BigTracksMC>(), candLambdaC.prong2_as<aod::BigTracksMC>()};
                int8_t sign = 0;
                int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
                auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
                // Get the corresponding MC particle for soft pion
                auto particleSoftPi = candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle();

                //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
                //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.);
                double ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                double etaSigmaC(candSigmaC.eta()), etaLambdaC(candLambdaC.eta());
                double phiSigmaC(candSigmaC.phi()), phiLambdaC(candLambdaC.phi());
                double ptSoftPi(candSigmaC.prong1().pt()), etaSoftPi(candSigmaC.prong1().eta()), phiSoftPi(candSigmaC.prong1().phi());
                double ptGenSigmaC(particleSigmaC.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
                int origin = candSigmaC.originMcRec();
                auto channel = candLambdaC.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
                
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kProton ) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for reconstructed Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++ signal

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kPiPlus ) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for reconstructed Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0,++ signal

                } /// end candidate Λc+ → π+K-p (and charge conjugate)
            } /// end reconstructed Σc0 signal

            /// Reconstructed Σc++ signal
            else if(std::abs(candSigmaC.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi && (std::abs(chargeSigmaC) == 2)) {
                // Get the corresponding MC particle for SigmaC, found as the mother of the soft pion
                auto indexMcSigmaCRec = RecoDecay::getMother(particlesMc, candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen>>(), pdg::Code::kSigmaC0, true);
                auto particleSigmaC = particlesMc.rawIteratorAt(indexMcSigmaCRec);
                // Get the corresponding MC particle for Lc
                auto arrayDaughtersLc = array{candLambdaC.prong0_as<aod::BigTracksMC>(), candLambdaC.prong1_as<aod::BigTracksMC>(), candLambdaC.prong2_as<aod::BigTracksMC>()};
                int8_t sign = 0;
                int indexMcLcRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
                auto particleLambdaC = particlesMc.rawIteratorAt(indexMcLcRec);
                // Get the corresponding MC particle for soft pion
                auto particleSoftPi = candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle();

                //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
                //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.);
                double ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                double etaSigmaC(candSigmaC.eta()), etaLambdaC(candLambdaC.eta());
                double phiSigmaC(candSigmaC.phi()), phiLambdaC(candLambdaC.phi());
                double ptSoftPi(candSigmaC.prong1().pt()), etaSoftPi(candSigmaC.prong1().eta()), phiSoftPi(candSigmaC.prong1().phi());
                double ptGenSigmaC(particleSigmaC.pt()), ptGenLambdaC(particleLambdaC.pt()), ptGenSoftPi(particleSoftPi.pt());
                int origin = candSigmaC.originMcRec();
                auto channel = candLambdaC.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±

                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kProton ) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for reconstructed Σc++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0,++ signal

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && std::abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kPiPlus) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for reconstructed Σc++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroPlusPlusSig"), ptSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSigmaCZeroPlusPlusSig"), ptGenSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroPlusPlusSig"), etaSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroPlusPlusSig"), phiSigmaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSig"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptSigmaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptSigmaC, (double) channel);   // Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroPlusPlusSig"), ptSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenSoftPiSigmaCZeroPlusPlusSig"), ptGenSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroPlusPlusSig"), etaSoftPi, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroPlusPlusSig"), phiSoftPi, (double)origin, (double) channel);   // π ← Σc0,++ signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroPlusPlusSig"), ptLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPtGenLambdaCFromSigmaCZeroPlusPlusSig"), ptGenLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroPlusPlusSig"), etaLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroPlusPlusSig"), phiLambdaC, (double)origin, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSig"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigPrompt"), deltaMass, ptLambdaC, (double) channel);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroPlusPlusSigNonPrompt"), deltaMass, ptLambdaC, (double) channel);   // Λc+ ← Σc0,++ signal

                } /// end candidate Λc+ → π+K-p (and charge conjugate)
            } /// end reconstructed Σc++ signal

        } /// end loop on reconstructed Σc0,++

    }; /// end processMc
    PROCESS_SWITCH(HfTaskSigmaC, processMc, "Process MC", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSigmaC>(cfgc, TaskName{"hf-task-sigmac"})};
}