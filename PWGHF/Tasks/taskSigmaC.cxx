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
        {"RecoData/hPtSigmaCZero", "#Sigma_{c}^{0} candidates; #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaSigmaCZero", "#Sigma_{c}^{0} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiSigmaCZero", "#Sigma_{c}^{0} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassSigmaCZero", "#Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Σc++
        {"RecoData/hPtSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassSigmaCPlusPlus", "#Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Σc0,++
        {"RecoData/hPtSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassSigmaCZeroPlusPlus", "#Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0
        {"RecoData/hPtLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCZero", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc++
        {"RecoData/hPtLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{p}_{T}(#Sigma_{c}^{++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}},
        /// Λc+ ← Σc0,++
        {"RecoData/hPtLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{p}_{T}(#Sigma_{c}^{0,++}) (GeV/#it{c}); entries;", {HistType::kTH1D, {{36, 0., 36.}}}},
        {"RecoData/hEtaLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #eta; entries;", {HistType::kTH1D, {{36, -2., 2.}}}},
        {"RecoData/hPhiLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #varphi; entries;", {HistType::kTH1D, {{72, 0, 2*M_PI}}}},
        {"RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus", "#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++} candidates; #it{M}(pK#pi#pi) - #it{M}(pK#pi) (GeV/#it{c}^{2}); #it{p}_{T}(#Lambda_{c}^{+} #leftarrow #Sigma_{c}^{0,++}) (GeV/#it{c});", {HistType::kTH2F, {{200, 0.13, 0.23}, {36, 0., 36.}}}}}
    };

    /// @brief init function, to define the additional analysis histograms
    /// @param  
    void init(InitContext&) {
        /// TO DO: add histograms for MC, when required
        /// [...]
        if(doprocessMc) {
            
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
            /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
            if(isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
                massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                massLambdaC = invMassLcToPKPi(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("RecoData/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZero"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZero"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZero"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZero"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("RecoData/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                }
            } /// end candidate Λc+ → pK-π+ (and charge conjugate)
            /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
            if(isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
                massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                massLambdaC = invMassLcToPiKP(candLambdaC);
                deltaMass = massSigmaC - massLambdaC;
                /// fill the histograms
                if(chargeSigmaC == 0) {
                    registry.fill(HIST("RecoData/hPtSigmaCZero"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZero"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZero"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZero"), deltaMass, ptSigmaC);   // Σc0
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZero"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZero"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZero"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZero"), deltaMass, ptLambdaC);   // Λc+ ← Σc0
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
                } else { /// candidates with charge ++ (or --). Possible unexpected candidates with charge + (or -) already discared in HFCandidateCreatorSigmaCZeroPlusPlus.cxx
                    registry.fill(HIST("RecoData/hPtSigmaCPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCPlusPlus"), deltaMass, ptSigmaC);   // Σc++
                    registry.fill(HIST("RecoData/hPtSigmaCZeroPlusPlus"), ptSigmaC);
                    registry.fill(HIST("RecoData/hEtaSigmaCZeroPlusPlus"), etaSigmaC);
                    registry.fill(HIST("RecoData/hPhiSigmaCZeroPlusPlus"), phiSigmaC);
                    registry.fill(HIST("RecoData/hDeltaMassSigmaCZeroPlusPlus"), deltaMass, ptSigmaC);   // Σc0,++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc++
                    registry.fill(HIST("RecoData/hPtLambdaCFromSigmaCZeroPlusPlus"), ptLambdaC);
                    registry.fill(HIST("RecoData/hEtaLambdaCFromSigmaCZeroPlusPlus"), etaLambdaC);
                    registry.fill(HIST("RecoData/hPhiLambdaCFromSigmaCZeroPlusPlus"), phiLambdaC);
                    registry.fill(HIST("RecoData/hDeltaMassLambdaCFromSigmaCZeroPlusPlus"), deltaMass, ptLambdaC);   // Λc+ ← Σc0,++
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
            const auto& candLambdaC = candSigmaC.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
            const int isCandPKPiPiKP = isDecayToPKPiToPiKP(candLambdaC, candSigmaC);

            /// Reconstructed Σc0 signal
            if(std::abs(candSigmaC.flagMcMatchRec()) == 1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi && (chargeSigmaC == 0)) {
                // Get the corresponding MC particle, found as the mother of the soft pion
                auto indexMother = RecoDecay::getMother(particlesMc, candSigmaC.prong1_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandSigmaCMcGen>>(), pdg::Code::kSigmaC0, true);
                auto particleMother = particlesMc.rawIteratorAt(indexMother);
                //registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT
                
                //const int isCandLambdaCpKpi = (candLambdaC.isSelLcToPKPi() >= 1) && candSigmaC.statusSpreadLcMinvPKPiFromPDG(); // Λc+ → pK-π+ and within the requested mass to build the Σc0,++
                //const int isCandLambdaCpiKp = (candLambdaC.isSelLcToPiKP() >= 1) && candSigmaC.statusSpreadLcMinvPiKPFromPDG(); // Λc+ → π+K-p and within the requested mass to build the Σc0,++
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.), ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                auto ptSc = candSigmaC.pt();
                auto ptLc = candLambdaC.pt();
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if(isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for Σc0 and Σc0,++
                    /// TODO [...]

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if(isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for Σc0 and Σc0,++
                    /// TODO [...]

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
                double massSigmaC(-1.), massLambdaC(-1.), deltaMass(-1.), ptSigmaC(candSigmaC.pt()), ptLambdaC(candLambdaC.pt());
                auto ptSc = candSigmaC.pt();
                auto ptLc = candLambdaC.pt();
                /// candidate Λc+ → pK-π+ (and charge conjugate) within the range of M(pK-π+) chosen in the Σc0,++ builder
                if(isCandPKPiPiKP == 1 || isCandPKPiPiKP == 3) {
                    massSigmaC = invMassScRecoLcToPKPi(candSigmaC);
                    massLambdaC = invMassLcToPKPi(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;

                    /// Fill the histograms for Σc0 and Σc0,++
                    /// TODO [...]

                } /// end candidate Λc+ → pK-π+ (and charge conjugate)
                /// candidate Λc+ → π+K-p (and charge conjugate) within the range of M(π+K-p) chosen in the Σc0,++ builder
                if(isCandPKPiPiKP == 2 || isCandPKPiPiKP == 3) {
                    massSigmaC = invMassScRecoLcToPiKP(candSigmaC);
                    massLambdaC = invMassLcToPiKP(candLambdaC);
                    deltaMass = massSigmaC - massLambdaC;
                    
                    /// Fill the histograms for Σc0 and Σc0,++
                    /// TODO [...]

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