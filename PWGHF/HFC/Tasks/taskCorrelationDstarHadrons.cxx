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

/// \file taskCorrelationrDstarHadron.cxx   
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Shyam Kumar <shyam.kumar@cern.ch>

// Framework
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// PWGHF
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"


using namespace o2;
// using namespace o2::dstar;
using namespace o2::framework;
using namespace o2::framework::expressions;

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringPoolBin = "Pool Bin Number;";



namespace o2::dstar::correlation{
    const int nBinsPtCorrelation = 8;

    const double binsPtCorrelations[nBinsPtCorrelation + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
    auto vecBinsPtCorrelations = std::vector<double>{binsPtCorrelations, binsPtCorrelations + nBinsPtCorrelation + 1};

    const double signalRegionLefBoundDefault[nBinsPtCorrelation] = {1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490};
    auto vecSignalRegionLefBoundDefault = std::vector<double>{signalRegionLefBoundDefault , signalRegionLefBoundDefault + nBinsPtCorrelation };

    const double signalRegionRightBoundDefault[nBinsPtCorrelation] = {1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890};
    auto vecSignalRegionRightBoundDefault = std::vector<double>{signalRegionRightBoundDefault, signalRegionRightBoundDefault + nBinsPtCorrelation };

    const double sidebandLeftOuterDefault[nBinsPtCorrelation] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
    auto vecSidebandLeftOuterDefault = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nBinsPtCorrelation }; 

    const double sidebandLeftInnerDefault[nBinsPtCorrelation] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
    auto vecSidebandLeftInnerDefault = std::vector<double>{sidebandLeftInnerDefault , sidebandLeftInnerDefault + nBinsPtCorrelation };

    const double sidebandRightInnerDefault[nBinsPtCorrelation] = {1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130};
    auto vecSidebandRightInnerDefault = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nBinsPtCorrelation };

    const double sidebandRightOuterDefault[nBinsPtCorrelation] = {1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690};
    auto vecSidebandRightOuterDefault = std::vector<double>{sidebandRightOuterDefault , sidebandRightOuterDefault + nBinsPtCorrelation };

    const int npTBinsEfficiency = o2::analysis::hf_cuts_dstar_to_d0_pi::nBinsPt;
    std::vector<double> efficiencyDstar(npTBinsEfficiency); // line # 76 in taskCorrelationDstarHadron.cxx; why (npTBinsEfficiency+1) ?
}

using namespace o2::dstar;

// Dstar-Hadron correlation pair
struct HfTaskCorrelationDstarHadron{

    Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
    // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_to_pi_k_pi (i.e. the mass pT bins), but can be redefined via json files
    Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{correlation::vecBinsPtCorrelations}, "pT bin limits for correlation plots"};
    Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits for efficiency"};
    Configurable<std::vector<double>> efficiencyDstar{"efficiencyDstar", std::vector<double>{correlation::efficiencyDstar}, "efficiency values for Dstar vs pT bin"};

    Configurable<std::vector<double>> signalRegionLefBound{"signalRegionLefBound",std::vector<double>{correlation::vecSignalRegionLefBoundDefault},"left boundary of signal region vs pT"};
    Configurable<std::vector<double>> signalRegionRightBound{"signalRegionRightBound",std::vector<double>{correlation::vecSignalRegionRightBoundDefault},"right boundary of signal region vs pT"};
    Configurable<std::vector<double>> leftSidebandOuterBoundary{"leftSidebandOuterBoundary",std::vector<double>{correlation::vecSidebandLeftOuterDefault},"left sideband outer boundary vs pT"};
    Configurable<std::vector<double>> leftSidebandInnerBoundary{"leftSidebandInnerBoundary", std::vector<double>{correlation::vecSidebandLeftInnerDefault},"left sideband inner boundary vs pT"};
    Configurable<std::vector<double>> rightSidebandOuterBoundary{"rightSidebandOuterBoundary",std::vector<double>{correlation::vecSidebandRightOuterDefault},"right sideband outer baoundary vs pT"};
    Configurable<std::vector<double>> rightSidebandInnerBoundary{"rightSidebandInnerBoundary",std::vector<double>{correlation::vecSidebandRightInnerDefault},"right sideband inner boundary"};



    HistogramRegistry registry{
        "registry",
        {
            // {"hCorrel2DVsPtSignalRegion",stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron +stringPoolBin + "entries",{HistType::kTHnSparseD,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2},{50,0.,50.},{11,0.,11.},{9,0.,9.}}}},
            // {"hCorrel2DPtIntSignalRegion",stringDHadron + stringSignal + stringDeltaPhi +stringDeltaEta + "entries",{HistType::kTH2D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2}}}},
            // {"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries",{HistType::kTH1D,{{40,-2.,2}}}},
            // {"hDeltaPhiPtIntSignalRegion",stringDHadron + stringSignal + stringDeltaPhi + "entries",{HistType::kTH1D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf}}}},
            // {"hCorrel2DVsPtSidebands",stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries",{HistType::kTHnSparseD,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2},{50,0.,50.},{11,0.,11.},{9,0.,9.}}}},
            // {"hCorrel2DPtIntSidebands",stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries",{HistType::kTH2D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2}}}},
            // {"hDeltaEtaPtIntSidebands",stringDHadron + stringSideband + stringDeltaEta + "entries",{HistType::kTH1D,{{40,-2.,2}}}},
            // {"hDeltaPhiPtIntSidebands",stringDHadron + stringSideband + stringDeltaPhi + "entries",{HistType::kTH1D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf}}}}
        }

    };
    
    void init(InitContext &){

        auto axisPtDstar = (std::vector<double>)binsPtEfficiency;
        AxisSpec axisSpecPtDstar = {axisPtDstar};

        registry.add("hCorrel2DVsPtSignalRegion",stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron +stringPoolBin + "entries",{HistType::kTHnSparseD,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2},axisSpecPtDstar,{11,0.,11.},{9,0.,9.}}}, true);
        registry.add("hCorrel2DPtIntSignalRegion",stringDHadron + stringSignal + stringDeltaPhi +stringDeltaEta + "entries",{HistType::kTH2D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2}}}, true);
        registry.add("hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries",{HistType::kTH1D,{{40,-2.,2}}}, true);
        registry.add("hDeltaPhiPtIntSignalRegion",stringDHadron + stringSignal + stringDeltaPhi + "entries",{HistType::kTH1D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf}}}, true);
        registry.add("hCorrel2DVsPtSidebands",stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries",{HistType::kTHnSparseD,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2},axisSpecPtDstar,{11,0.,11.},{9,0.,9.}}}, true);
        registry.add("hCorrel2DPtIntSidebands",stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries",{HistType::kTH2D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf},{40,-2.,2}}}, true);
        registry.add("hDeltaEtaPtIntSidebands",stringDHadron + stringSideband + stringDeltaEta + "entries",{HistType::kTH1D,{{40,-2.,2}}}, true);
        registry.add("hDeltaPhiPtIntSidebands",stringDHadron + stringSideband + stringDeltaPhi + "entries",{HistType::kTH1D,{{64,-o2::constants::math::PIHalf, 3.*o2::constants::math::PIHalf}}}, true);
    }

    void processData(aod::DstarHadronPair const & dstarHPairs){
        for(const auto & dstarHPair:dstarHPairs){
            float deltaPhi = dstarHPair.deltaPhi();
            float deltaEta = dstarHPair.deltaEta();
            float ptDstar = dstarHPair.ptDstar();
            float ptTrack = dstarHPair.ptTrack();
            int poolBin = dstarHPair.poolBin();
            float massDstar = dstarHPair.mDstar();

            int effBinPtDstar = o2::analysis::findBin(binsPtEfficiency,ptDstar);
            int corrBinPtDstar = o2::analysis::findBin(binsPtCorrelations,ptDstar);

            // reject candidate if outside pT ranges of interst
            if(corrBinPtDstar < 0 || effBinPtDstar < 0){
                continue;
            }
            if(ptTrack  > 10.0){
                ptTrack = 10.5;
            }
            float netEfficiencyWeight = 1.0;
            float efficiencyWeightTracks = 1.0;

            if(applyEfficiency){
                netEfficiencyWeight = 1.0 / (efficiencyDstar->at(effBinPtDstar) * efficiencyWeightTracks);
            }

             // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
            if(massDstar > signalRegionLefBound->at(ptDstar) && massDstar < signalRegionRightBound->at(ptDstar)){
                // in signal region
                registry.fill(HIST("hCorrel2DVsPtSignalRegion"),deltaPhi, deltaEta, ptDstar , ptTrack, poolBin, netEfficiencyWeight);
                registry.fill(HIST("hCorrel2DPtIntSignalRegion"),deltaPhi, deltaEta, netEfficiencyWeight);
                registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, netEfficiencyWeight);
                registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, netEfficiencyWeight);
            }else if((massDstar > leftSidebandOuterBoundary->at(ptDstar) && massDstar < leftSidebandInnerBoundary->at(ptDstar))
                ||(massDstar > rightSidebandInnerBoundary->at(ptDstar) && massDstar < rightSidebandOuterBoundary->at(ptDstar))){
                registry.fill(HIST("hCorrel2DVsPtSidebands"),deltaPhi, deltaEta, ptDstar , ptTrack, poolBin, netEfficiencyWeight);
                registry.fill(HIST("hCorrel2DPtIntSidebands"),deltaPhi, deltaEta, netEfficiencyWeight);
                registry.fill(HIST("hDeltaEtaPtIntSidebands"),deltaEta, netEfficiencyWeight);
                registry.fill(HIST("hDeltaPhiPtIntSidebands"),deltaPhi, netEfficiencyWeight);
            }

        }

    }
    PROCESS_SWITCH(HfTaskCorrelationDstarHadron,processData," process data only", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDstarHadron>(cfgc)};
}
