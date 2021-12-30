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
///
/// \brief Task to check PID efficiency
/// \author Łukasz Sawicki
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/MC.h"
#include <TParameter.h>


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PIDV2 {
    double combindeSignal(float val1, float val2)
    {
        return sqrt(pow(val1, 2) + pow(val2, 2));
    }

    int indexofSmallestElement(const float array[], int size)
    {
        int index = 0;
        for(int i = 1; i < size; i++)
        {
            if(array[i] < array[index])
                index = i;              
        }
        return index;
    }

    void pid_false(int i, int pdgCode, double pt) {
        // check contamination
        if (i == 0) {
            if (pdgCode == 211)
                contamination_registry.fill(HIST("211in211"), pt);
            if (pdgCode == -211)
                contamination_registry.fill(HIST("0211in211"), pt);

            if (pdgCode == 2212)
                contamination_registry.fill(HIST("2212in211"), pt);
            if (pdgCode == -2212)
                contamination_registry.fill(HIST("02212in211"), pt);

            if (pdgCode == 321)
                contamination_registry.fill(HIST("321in211"), pt);
            if (pdgCode == -321)
                contamination_registry.fill(HIST("0321in211"), pt);

            if (pdgCode == 11)
                contamination_registry.fill(HIST("11in211"), pt);
            if (pdgCode == -11)
                contamination_registry.fill(HIST("011in211"), pt);

            if (pdgCode == 13)
                contamination_registry.fill(HIST("13in211"), pt);
            if (pdgCode == -13)
                contamination_registry.fill(HIST("013in211"), pt);
        }
        else if (i == 1) {
            if (pdgCode == 211)
                contamination_registry.fill(HIST("211in2212"), pt);
            if (pdgCode == -211)
                contamination_registry.fill(HIST("0211in2212"), pt);

            if (pdgCode == 2212)
                contamination_registry.fill(HIST("2212in221"), pt);
            if (pdgCode == -2212)
                contamination_registry.fill(HIST("02212in221"), pt);

            if (pdgCode == 321)
                contamination_registry.fill(HIST("321in2212"), pt);
            if (pdgCode == -321)
                contamination_registry.fill(HIST("0321in2212"), pt);

            if (pdgCode == 11)
                contamination_registry.fill(HIST("11in2212"), pt);
            if (pdgCode == -11)
                contamination_registry.fill(HIST("011in2212"), pt);

            if (pdgCode == 13)
                contamination_registry.fill(HIST("13in2212"), pt);
            if (pdgCode == -13)
                contamination_registry.fill(HIST("013in2212"), pt);
        }
        else if (i == 2) {
            if (pdgCode == 211)
                contamination_registry.fill(HIST("211in321"), pt);
            if (pdgCode == -211)
                contamination_registry.fill(HIST("0211in321"), pt);

            if (pdgCode == 2212)
                contamination_registry.fill(HIST("2212in321"), pt);
            if (pdgCode == -2212)
                contamination_registry.fill(HIST("02212in321"), pt);

            if (pdgCode == 321)
                contamination_registry.fill(HIST("321in321"), pt);
            if (pdgCode == -321)
                contamination_registry.fill(HIST("0321in321"), pt);

            if (pdgCode == 11)
                contamination_registry.fill(HIST("11in321"), pt);
            if (pdgCode == -11)
                contamination_registry.fill(HIST("011in321"), pt);

            if (pdgCode == 13)
                contamination_registry.fill(HIST("13in321"), pt);
            if (pdgCode == -13)
                contamination_registry.fill(HIST("013in321"), pt);
        }
    }

    // Histograms names, always keep the same order!
    static const int numParticles = 10;

    static constexpr std::string_view pid_true_registryH[numParticles] = {"pid_true_211", "pid_true_2212", "pid_true_321", "pid_true_11", "pid_true_13", "pid_true_0211", "pid_true_02212", "pid_true_0321", "pid_true_011", "pid_true_013"};
    static constexpr std::string_view pid_false_registryH[numParticles] = {"pid_false_211", "pid_false_2212", "pid_false_321", "pid_false_11", "pid_false_13", "pid_false_0211", "pid_false_02212", "pid_false_0321", "pid_false_011", "pid_false_013"};

    static constexpr std::string_view tpc_signal_registry_trueH[numParticles] = {"tpc_signal_pid_true_211", "tpc_signal_pid_true_2212", "tpc_signal_pid_true_321", "tpc_signal_pid_true_11", "tpc_signal_pid_true_13", "tpc_signal_pid_true_0211", "tpc_signal_pid_true_02212", "tpc_signal_pid_true_0321", "tpc_signal_pid_true_011", "tpc_signal_pid_true_013"};
    static constexpr std::string_view tpc_signal_registry_falseH[numParticles] = {"tpc_signal_pid_false_211", "tpc_signal_pid_false_2212", "tpc_signal_pid_false_321", "tpc_signal_pid_false_11", "tpc_signal_pid_false_13", "tpc_signal_pid_false_0211", "tpc_signal_pid_false_02212", "tpc_signal_pid_false_0321", "tpc_signal_pid_false_011", "tpc_signal_pid_false_013"};
    
    static constexpr std::string_view tof_signal_registry_trueH[numParticles] = {"tof_signal_pid_true_211", "tof_signal_pid_true_2212", "tof_signal_pid_true_321", "tof_signal_pid_true_11", "tof_signal_pid_true_13", "tof_signal_pid_true_0211", "tof_signal_pid_true_02212", "tof_signal_pid_true_0321", "tof_signal_pid_true_011", "tof_signal_pid_true_013"};
    static constexpr std::string_view tof_signal_registry_falseH[numParticles] = {"tof_signal_pid_false_211", "tof_signal_pid_false_2212", "tof_signal_pid_false_321", "tof_signal_pid_false_11", "tof_signal_pid_false_13", "tof_signal_pid_false_0211", "tof_signal_pid_false_02212", "tof_signal_pid_false_0321", "tof_signal_pid_false_011", "tof_signal_pid_false_013"};

    static constexpr int pdgCodes[numParticles] = {211, 2212, 321, 11, 13, -211, -2212, -321, -11, -13};
    static constexpr int charges[numParticles] = {1, 1, 1, -1, -1, -1, -1, -1, 1, 1};
    static constexpr float ptCut[numParticles] = {0.5, 0.8, 0.5, 0.5, 0.5, 0.5, 0.8, 0.5, 0.5, 0.5};

    static const int max_p = 5;
    static const int nb_bins = 100;
    static const int nb_bins2D = 1000;

    // registry to keep particles from MC data
    HistogramRegistry mc_registry{
        "mc_registry",
        {
            {"mc_211", "MC Pions^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_11", "MC Electrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_321", "MC Kaons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_2212", "MC Protons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_13", "MC Muons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_0211", "MC Pions^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_011", "MC Positrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_0321", "MC Kaons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_02212", "MC Anitprotons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"mc_013", "MC Muons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"else", "MC else;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}}
        }
    };

    HistogramRegistry pid_true_registry{
        "pid_true_registry",
        {
            {"pid_true_211", "PID true Pions^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_11", "PID true Electrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_321", "PID true Kaons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_2212", "PID true Protons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_13", "PID true Muons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_0211", "PID true Pions^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_011", "PID true Positrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_0321", "PID true Kaons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_02212", "PID true Antiprotons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_true_013", "PID true Muons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}}
        }
    };

    HistogramRegistry pid_false_registry{
        "pid_false_registry",
        {
            {"pid_false_211", "PID false Pions^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_11", "PID false Electrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_321", "PID false Kaons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_2212", "PID false Protons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_13", "PID false Muons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_0211", "PID false Pions^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_011", "PID false Positrons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_0321", "PID false Kaons^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_02212", "PID false Antiprotons;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"pid_false_013", "PID false Muons^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}}
        }
    };

    HistogramRegistry contamination_registry{
        "contamination_registry",
        {
            {"211in211", "Pions^{+} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"2212in211", "Protons in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"321in211", "Kaons^{+} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"11in211", "Electrons in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"13in211", "Muons^{-} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},

            {"0211in211", "Pions^{-} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"02212in211", "Antiprotons in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"0321in211", "Kaons^{-} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"011in211", "Positrons in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"013in211", "Muons^{+} in Pions^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},

            {"211in2212", "Pions^{+} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"2212in2212", "Protons in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"321in2212", "Kaons^{+} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"11in2212", "Electrons in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"13in2212", "Muons^{-} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},

            {"0211in2212", "Pions^{-} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"02212in2212", "Antiprotons in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"0321in2212", "Kaons^{-} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"011in2212", "Positrons in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"013in2212", "Muons^{+} in Protons sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},

            {"211in321", "Pions^{+} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"2212in321", "Protons in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"321in321", "Kaons^{+} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"11in321", "Electrons in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"13in321", "Muons^{-} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},

            {"0211in321", "Pions^{-} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"02212in321", "Antiprotons in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"0321in321", "Kaons^{-} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"011in321", "Positrons in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}},
            {"013in321", "Muons^{+} in Kaons^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{nb_bins, 0, max_p}}}}
        }
    };

    HistogramRegistry tpc_signal_registry{
        "tpc_signal_registry",
        {
            {"tpc_signal_PID_true", "TPCSignal PID true;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_PID_false", "TPCSignal PID false;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},

            {"tpc_signal_pid_true_211", "TPCSignal PID true #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_11", "TPCSignal PID true e^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_321", "TPCSignal PID true K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_2212", " TPCSignal PID true p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_13", "TPCSignal PID true #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_211", "TPCSignal PID false #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_11", "TPCSignal PID false e^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_321", "TPCSignal PID false K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_2212", "TPCSignal PID false p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_13", "TPCSignal PID false #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},

            {"tpc_signal_pid_true_0211", "TPCSignal PID true #pi^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_011", "TPCSignal PID true e^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_0321", "TPCSignal PID true K^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_02212", " TPCSignal PID true #bar{p};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_true_013", "TPCSignal PID true #mu^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_0211", "TPCSignal PID false #pi^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_011", "TPCSignal PID false e^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_0321", "TPCSignal PID false K^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_02212", "TPCSignal PID false #bar{p};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}},
            {"tpc_signal_pid_false_013", "TPCSignal PID false #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{nb_bins2D, 0, 10}, {nb_bins2D, 0, 300}}}}
        }
    };

    HistogramRegistry tof_signal_registry {
        "tof_signal_registry",
        {
            {"tof_signal_PID_true", "TOF PID true;p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_PID_false", "TOF PID false;p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},

            {"tof_signal_pid_true_211", "TOF PID true #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_2212", "TOF PID true p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_321", "TOF PID true K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_13", "TOF PID true #mu^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_11", "TOF PID true e^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_211", "TOF PID false #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_2212", "TOF PID false p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_321", "TOF PID false K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_13", "TOF PID false #mu^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_11", "TOF PID false e^{-};p (GeV/c); TOF b#eta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},

            {"tof_signal_pid_true_0211", "TOF PID true #pi^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_02212", "TOF PID true #bar{p};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_0321", "TOF PID true K^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_013", "TOF PID true #mu^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_true_011", "TOF PID true e^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_0211", "TOF PID false #pi^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_02212", "TOF PID false #bar{p};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_0321", "TOF PID false K^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_013", "TOF PID false #mu^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}},
            {"tof_signal_pid_false_011", "TOF PID false e^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{nb_bins2D, 0.2, 10}, {110, 0, 1.1}}}}
        }
    };

    template <std::size_t i, typename T>
    void PID_simple(const T& track, const int pdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arr_len) {
        const float p = track.p();

        bool false_pid = false;

        if ((p < ptCut[i]) & (track.sign() == charges[i])) {
            if (abs(tpcNSigmas[i]) < nsigmacut.value) {
                if (pdgCode == pdgCodes[i]) {
                    pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                    tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                    tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                    tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                    tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                } else {
                    false_pid = true;
                }
            }
        } 
        else if ((p >= ptCut[i]) & (track.sign() == charges[i])) {
            if (sqrt(pow(tpcNSigmas[i], 2) + pow(tofNSigmas[i], 2)) < nsigmacut.value) {
                if (pdgCode == pdgCodes[i]) {
                    pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                    tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                    tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                    tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                    tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                } else {
                    false_pid = true;
                }
            }
        }

        if (false_pid) {
            pid_false_registry.fill(HIST(pid_false_registryH[i]), track.pt());
            tpc_signal_registry.fill(HIST(tpc_signal_registry_falseH[i]), track.p(), track.tpcSignal());
            tpc_signal_registry.fill(HIST("tpc_signal_PID_false"), track.p(), track.tpcSignal());
            tof_signal_registry.fill(HIST("tof_signal_PID_false"), track.p(), track.beta());
            tof_signal_registry.fill(HIST(tof_signal_registry_falseH[i]), track.p(), track.beta());

            double pt = track.pt();
            pid_false(i, pdgCode, pt);
        }
    }

    template <std::size_t i, typename T>
    void PID_min_strategy(const T& track, const int pdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arr_len) {
        const float p = track.p();

        // list of Nsigmas for particles
        float particleNSigma[arr_len] = {};

        bool false_pid = false;

        // calculate Nsigmas for every particle
        for (int j = 0; j < arr_len; ++j) {
            if (p < ptCut[j]) {
                particleNSigma[j] = abs(tpcNSigmas[j]);
            } else if (p >= ptCut[j]) {
                particleNSigma[j] = combindeSignal(tpcNSigmas[j], tofNSigmas[j]);
            }
        }

        if ((p < ptCut[i]) & (track.sign() == charges[i])) {
            float tmp_NSigma = abs(tpcNSigmas[i]);
            if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arr_len) == i)) {
                if (pdgCode == pdgCodes[i]) {
                    pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                    tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                    tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                    tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                    tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                } else {
                    false_pid = true;
                }
            }
        }
        else if ((p >= ptCut[i]) & (track.sign() == charges[i])) {
            float tmp_NSigma = combindeSignal(tpcNSigmas[i], tofNSigmas[i]);
            if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arr_len) == i)) {
                if (pdgCode == pdgCodes[i]) {
                    pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                    tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                    tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                    tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                    tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                } else {
                    false_pid = true;
                }
            }
        }

        if (false_pid) {
            pid_false_registry.fill(HIST(pid_false_registryH[i]), track.pt());
            tpc_signal_registry.fill(HIST(tpc_signal_registry_falseH[i]), track.p(), track.tpcSignal());
            tpc_signal_registry.fill(HIST("tpc_signal_PID_false"), track.p(), track.tpcSignal());
            tof_signal_registry.fill(HIST("tof_signal_PID_false"), track.p(), track.beta());
            tof_signal_registry.fill(HIST(tof_signal_registry_falseH[i]), track.p(), track.beta());

            double pt = track.pt();
            pid_false(i, pdgCode, pt);
        }
    }

    template <std::size_t i, typename T>
    void PID_exclusive_strategy(const T& track, const int pdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arr_len) {
        const float p = track.p();
        // const int pdgCode = track.pdgCode();

        bool false_pid = false;

        // list of Nsigmas for particles
        float particleNSigma[arr_len] = {};

        // calculate Nsigmas for every particle
        for (int j = 0; j < arr_len; ++j) {
            if (p < ptCut[j]) {
                particleNSigma[j] = abs(tpcNSigmas[j]);
            } else if (p >= ptCut[j]) {
                particleNSigma[j] = combindeSignal(tpcNSigmas[j], tofNSigmas[j]);
            }
        }

        // check how many particles satisfy the condition
        int counts = 0;
        for (int j = 0; j < arr_len; ++j) {
            if (particleNSigma[j] < nsigmacut.value) {
                counts ++;
            }
        }

        if (counts == 1) {
            if ((p < ptCut[i]) & (track.sign() == charges[i])) {
                float tmp_NSigma = abs(tpcNSigmas[i]);
                if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arr_len) == i)) {
                    if (pdgCode == pdgCodes[i]) {
                        pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                        tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                        tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                        tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                        tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                    } else {
                        false_pid = true;
                    }
                }
            }
            else if ((p >= ptCut[i]) & (track.sign() == charges[i])) {
                float tmp_NSigma = combindeSignal(tpcNSigmas[i], tofNSigmas[i]);
                if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arr_len) == i)) {
                    if (pdgCode == pdgCodes[i]) {
                        pid_true_registry.fill(HIST(pid_true_registryH[i]), track.pt());
                        tpc_signal_registry.fill(HIST(tpc_signal_registry_trueH[i]), track.p(), track.tpcSignal());
                        tpc_signal_registry.fill(HIST("tpc_signal_PID_true"), track.p(), track.tpcSignal());
                        tof_signal_registry.fill(HIST("tof_signal_PID_true"), track.p(), track.beta());
                        tof_signal_registry.fill(HIST(tof_signal_registry_trueH[i]), track.p(), track.beta());
                    } else {
                        false_pid = true;
                    }
                }
            }
            if (false_pid) {
                pid_false_registry.fill(HIST(pid_false_registryH[i]), track.pt());
                tpc_signal_registry.fill(HIST(tpc_signal_registry_falseH[i]), track.p(), track.tpcSignal());
                tpc_signal_registry.fill(HIST("tpc_signal_PID_false"), track.p(), track.tpcSignal());
                tof_signal_registry.fill(HIST("tof_signal_PID_false"), track.p(), track.beta());
                tof_signal_registry.fill(HIST(tof_signal_registry_falseH[i]), track.p(), track.beta());
                
                double pt = track.pt();
                pid_false(i, pdgCode, pt);
            }
        } 
    }

    void init(o2::framework::InitContext&) {
        mc_registry.get<TH1>(HIST("mc_211")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_2212")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_321")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_11")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_13")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_0211")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_02212")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_0321")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_011")) -> Sumw2();
        mc_registry.get<TH1>(HIST("mc_013")) -> Sumw2();

        pid_true_registry.get<TH1>(HIST("pid_true_211")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_2212")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_321")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_11")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_13")) -> Sumw2();

        pid_true_registry.get<TH1>(HIST("pid_true_0211")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_02212")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_0321")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_011")) -> Sumw2();
        pid_true_registry.get<TH1>(HIST("pid_true_013")) -> Sumw2();

        pid_false_registry.get<TH1>(HIST("pid_false_211")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_2212")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_321")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_11")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_13")) -> Sumw2();

        pid_false_registry.get<TH1>(HIST("pid_false_0211")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_02212")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_0321")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_011")) -> Sumw2();
        pid_false_registry.get<TH1>(HIST("pid_false_013")) -> Sumw2();
    }

    Configurable<float> nsigmacut{"nsigmacut", 1.5, "Value of the NSigma cut"};
    Configurable<int> strategy{"Strategy", 3, "1-Simple, 2-Minimum, 3-Exclusive"};


    Filter trackFilter = aod::track::isGlobalTrack == static_cast<uint8_t>(true);
    using pidTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TracksExtended, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTPCEl, aod::pidTPCMu, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa, aod::pidTOFEl, aod::pidTOFMu>>;
    void process(pidTracks const& tracks, aod::McParticles const& mcParticles) {
        for (auto& track : tracks) {
            auto particle = track.mcParticle();
            int pdgCode = particle.pdgCode();;

            // pions
            if (pdgCode == 211) {
                mc_registry.fill(HIST("mc_211"), track.pt());
            } 
            else if (pdgCode == -211) {
                mc_registry.fill(HIST("mc_0211"), track.pt());
            } 
            // protons
            else if (pdgCode == 2212) {
                mc_registry.fill(HIST("mc_2212"), track.pt());
            }
            else if (pdgCode == -2212) {
                mc_registry.fill(HIST("mc_02212"), track.pt());
            } 
            //kaons
            else if (pdgCode == 321) {
                mc_registry.fill(HIST("mc_321"), track.pt());
            }
            else if (pdgCode == -321) {
                mc_registry.fill(HIST("mc_0321"), track.pt());
            }
            // electrons
            else if (pdgCode == 11) {
                mc_registry.fill(HIST("mc_11"), track.pt());;
            }
            else if (pdgCode == -11) {
                mc_registry.fill(HIST("mc_011"), track.pt());
            }
            // muons
            else if (pdgCode == 13) {
                mc_registry.fill(HIST("mc_13"), track.pt());
            }
            else if (pdgCode == -13) {
                mc_registry.fill(HIST("mc_013"), track.pt());
            }
            else {
                mc_registry.fill(HIST("else"), track.pt());
            }

            // PID
            // indexes for particles: 0-Pi, 1-Pr, 2-Ka, 3-El, 4-Mu
            const float tpcNSigmas[numParticles] = {track.tpcNSigmaPi(), track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaEl(), track.tpcNSigmaMu()};
            const float tofNSigmas[numParticles] = {track.tofNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaEl(), track.tofNSigmaMu()};

            if (strategy.value == 1) {
                // a) strategia 1. tylko cutNSigma
                PID_simple<0>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<1>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<2>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<3>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<4>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<5>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<6>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<7>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<8>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_simple<9>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
            } else if (strategy.value == 2) {
                // b) strategia 2. NCutSigma + warunek na minimalne NSigma
                PID_min_strategy<0>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<1>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<2>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<3>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<4>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<5>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<6>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<7>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<8>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_min_strategy<9>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
            } else if (strategy.value == 3) {
                // c) strategia 3. exclusive NSigma
                PID_exclusive_strategy<0>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<1>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<2>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<3>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<4>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<5>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<6>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<7>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<8>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
                PID_exclusive_strategy<9>(track, pdgCode, tpcNSigmas, tofNSigmas, 5);
            }
        } // for
    }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec {
        adaptAnalysisTask<PIDV2>(cfgc),
    };
}