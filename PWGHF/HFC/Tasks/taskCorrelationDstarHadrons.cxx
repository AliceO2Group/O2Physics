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

/// \file correlatorDstarHadron.cxx 
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


// Dstar-Hadron correlation pair
struct HfTaskCorrelationDplusHadrons{

    // string definitions, used for histogram axis labels
    const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
    const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
    const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
    const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
    const TString stringDHadron = "D,Hadron candidates ";
    const TString stringSignal = "signal region;";
    const TString stringSideband = "sidebands;";

    HistogramRegistry registry{
        "registry",
        {
            {"hCorrel2DVsPtSignalRegion",stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries",{HistType::kTHnSparseD,{}}},
            {"hCorrel2DPtIntSignalRegion"," ",{HistType::kTH2D,{}}},
            {"hDeltaEtaPtIntSignalRegion"," ",{HistType::kTH1D,{}}},
            {"hDeltaPhiPtIntSignalRegion"," ",{HistType::kTH1D,{}}},
            {"hCorrel2DVsPtSidebands"," ",{HistType::kTHnSparseD,{}}},
            {"hCorrel2DPtIntSidebands"," ",{HistType::kTH2D,{}}},
            {"hDeltaEtaPtIntSidebands"," ",{HistType::kTH1D,{}}},
            {"hDeltaPhiPtIntSidebands"," ",{HistType::kTH1D,{}}}
        }

    };

};
