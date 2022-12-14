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
        {"Data/hEtaSigmaCZero", "#Sigma_{c}^{0} candidates; #eta(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiSigmaCZero", "#Sigma_{c}^{0} candidates; #varphi(#Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCZero", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Σc++
        {"Data/hPtSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #eta(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #varphi(#Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Σc0,++
        {"Data/hPtSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #eta(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi(#Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        {"Data/hPtSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}}, /// soft π
        {"Data/hEtaSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}}, /// soft π
        {"Data/hPhiSoftPiSigmaCZeroPlusPlus", "#pi^{#pm} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}}, /// soft π
        /// Λc+ ← Σc0
        {"Data/hPtLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc++
        {"Data/hPtLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0,++
        {"Data/hPtLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"Data/hEtaLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"Data/hPhiLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}); entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"Data/hDeltaMassLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}
    };

    /// @brief init function, to define the additional analysis histograms
    /// @param  
    void init(InitContext&) {
        /// TO DO: add histograms for MC, when required
        /// [...]
        if(doprocessMc) {
            /// Reconstructed Σc0 signal
            registry.add("MC/reconstructed/hPtSigmaCZeroSig", "#Sigma_{c}^{0} signal; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hEtaSigmaCZeroSig", "#Sigma_{c}^{0} signal; #eta(#Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hPhiSigmaCZeroSig", "#Sigma_{c}^{0} signal; #varphi(#Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSig", "#Sigma_{c}^{0} signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt", "#Sigma_{c}^{0} prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt", "#Sigma_{c}^{0} non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hPtSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}}); /// soft π
            registry.add("MC/reconstructed/hEtaSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}}); /// soft π
            registry.add("MC/reconstructed/hPhiSoftPiSigmaCZeroSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{0} signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}}); /// soft π
            /// Reconstructed Σc++ signal
            registry.add("MC/reconstructed/hPtSigmaCPlusPlusSig", "#Sigma_{c}^{++} signal; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hEtaSigmaCPlusPlusSig", "#Sigma_{c}^{++} signal; #eta(#Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hPhiSigmaCPlusPlusSig", "#Sigma_{c}^{++} signal; #varphi(#Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig", "#Sigma_{c}^{++} signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt", "#Sigma_{c}^{++} prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt", "#Sigma_{c}^{++} non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} signal; #it{p}_{T}(#pi^{#pm} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}}); /// soft π
            registry.add("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} signal; #eta(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}}); /// soft π
            registry.add("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig", "#pi^{#pm} #leftarrow #Sigma_{c}^{++} signal; #varphi(#pi^{#pm} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}}); /// soft π
        
            /// Reconstructed Λc+ ← Σc0 signal
            registry.add("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            /// Reconstructed Λc+ ← Σc++ signal
            registry.add("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} signal; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, 0., 36.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} signal; #eta(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{36, -2., 2.}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} signal; #varphi(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}); origin (1 : prompt; 2: non-prompt);", {HistType::kTH2D, {{72, 0, 2*M_PI}, {2, 0.5, 2.5}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
            registry.add("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} non-prompt signal; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2D, {{200, 0.13, 0.23}, {36, 0., 36.}}});
        }
    }; /// end init

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

    /// @brief process function to fill the histograms needed in analysis (MC)
    /// @param candidatesSigmaC are the reconstructed candidate Σc0,++ with MC info
    /// @param particlesMc are the generated particles with flags wheter they are Σc0,++ or not
    /// @param 
    void processMc(const soa::Join<aod::HfCandSigmaC, aod::HfCandSigmaCMcRec>& candidatesSigmaC,
    soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen> const& particlesMc,
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

            /// Fill histograms
            /// TODO [...]

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

            /// Reconstructed Σc0 signal
            if(std::abs(candSigmaC.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi && (chargeSigmaC == 0)) {
                // Get the corresponding MC particle, found as the mother of the soft pion
                auto indexMother = RecoDecay::getMother(particlesMc, candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen>>(), pdg::Code::kSigmaC0, true);
                auto particleMother = particlesMc.rawIteratorAt(indexMother);
                //registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT
                
                //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
                //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.);
                double ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                double etaSigmaC(candSigmaC.eta()), etaLambdaC(candLambdaC.eta());
                double phiSigmaC(candSigmaC.phi()), phiLambdaC(candLambdaC.phi());
                double ptSoftPi(candSigmaC.prong1().pt()), etaSoftPi(candSigmaC.prong1().eta()), phiSoftPi(candSigmaC.prong1().phi());
                int origin = candSigmaC.originMcRec();
                
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kProton ) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for reconstructed Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroSig"), ptSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroSig"), etaSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroSig"), phiSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroSig"), ptSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroSig"), etaSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroSig"), phiSoftPi, origin);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig"), ptLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig"), etaLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig"), phiLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    /// [...]

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kPiPlus ) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for reconstructed Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCZeroSig"), ptSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCZeroSig"), etaSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCZeroSig"), phiSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCZeroSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCZeroSig"), ptSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCZeroSig"), etaSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCZeroSig"), phiSoftPi, origin);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCZeroSig"), ptLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCZeroSig"), etaLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCZeroSig"), phiLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCZeroSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    /// [...]

                } /// end candidate Λc+ → π+K-p (and charge conjugate)
            } /// end reconstructed Σc0 signal

            /// Reconstructed Σc++ signal
            else if(std::abs(candSigmaC.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi && (std::abs(chargeSigmaC) == 2)) {
                // Get the corresponding MC particle, found as the mother of the soft pion
                auto indexMother = RecoDecay::getMother(particlesMc, candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen>>(), pdg::Code::kSigmaCPlusPlus, true);
                auto particleMother = particlesMc.rawIteratorAt(indexMother);
                //registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT

                //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
                //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.);
                double ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                double etaSigmaC(candSigmaC.eta()), etaLambdaC(candLambdaC.eta());
                double phiSigmaC(candSigmaC.phi()), phiLambdaC(candLambdaC.phi());
                double ptSoftPi(candSigmaC.prong1().pt()), etaSoftPi(candSigmaC.prong1().eta()), phiSoftPi(candSigmaC.prong1().phi());
                int origin = candSigmaC.originMcRec();

                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) && abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kProton ) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for reconstructed Σc++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCPlusPlusSig"), ptSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCPlusPlusSig"), etaSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCPlusPlusSig"), phiSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig"), ptSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig"), etaSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig"), phiSoftPi, origin);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig"), ptLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig"), etaLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig"), phiLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    /// [...]

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if((isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) && abs(candLambdaC.prong0_as<aod::BigTracksMC>().mcParticle().pdgCode())==kPiPlus) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for reconstructed Σc++ signal
                    registry.fill(HIST("MC/reconstructed/hPtSigmaCPlusPlusSig"), ptSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSigmaCPlusPlusSig"), etaSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSigmaCPlusPlusSig"), phiSigmaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSig"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigPrompt"), deltaMass, ptSigmaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassSigmaCPlusPlusSigNonPrompt"), deltaMass, ptSigmaC);   // Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtSoftPiSigmaCPlusPlusSig"), ptSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaSoftPiSigmaCPlusPlusSig"), etaSoftPi, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiSoftPiSigmaCPlusPlusSig"), phiSoftPi, origin);   // π ← Σc0 signal
                    registry.fill(HIST("MC/reconstructed/hPtLambdaCFromSigmaCPlusPlusSig"), ptLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hEtaLambdaCFromSigmaCPlusPlusSig"), etaLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hPhiLambdaCFromSigmaCPlusPlusSig"), phiLambdaC, origin);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSig"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigPrompt"), deltaMass, ptLambdaC);
                    registry.fill(HIST("MC/reconstructed/hDeltaMassLambdaCFromSigmaCPlusPlusSigNonPrompt"), deltaMass, ptLambdaC);   // Λc+ ← Σc0 signal

                    /// Fill the histograms for reconstructed Σc0,++ signal
                    /// [...]

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