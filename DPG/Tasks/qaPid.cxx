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

struct pid {
  double combindeSignal(float val1, float val2)
  {
    return sqrt(pow(val1, 2) + pow(val2, 2));
  }

  int indexofSmallestElement(const float array[], int size)
  {
    int index = 0;
    for (int i = 1; i < size; i++) {
      if (array[i] < array[index])
        index = i;
    }
    return index;
  }

  void fillContaminationRegistry(int i, int PdgCode, double pt)
  {
    if (i == 0) {
      if (PdgCode == 211) {
        histReg.fill(HIST("contamination/211in211"), pt);
      }
      if (PdgCode == -211) {
        histReg.fill(HIST("contamination/0211in211"), pt);
      }

      if (PdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in211"), pt);
      }
      if (PdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in211"), pt);
      }

      if (PdgCode == 321) {
        histReg.fill(HIST("contamination/321in211"), pt);
      }
      if (PdgCode == -321) {
        histReg.fill(HIST("contamination/0321in211"), pt);
      }

      if (PdgCode == 11) {
        histReg.fill(HIST("contamination/11in211"), pt);
      }
      if (PdgCode == -11) {
        histReg.fill(HIST("contamination/011in211"), pt);
      }

      if (PdgCode == 13) {
        histReg.fill(HIST("contamination/13in211"), pt);
      }
      if (PdgCode == -13) {
        histReg.fill(HIST("contamination/013in211"), pt);
      }

    } else if (i == 1) {
      if (PdgCode == 211) {
        histReg.fill(HIST("contamination/211in2212"), pt);
      }
      if (PdgCode == -211) {
        histReg.fill(HIST("contamination/0211in2212"), pt);
      }

      if (PdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in221"), pt);
      }
      if (PdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in221"), pt);
      }

      if (PdgCode == 321) {
        histReg.fill(HIST("contamination/321in2212"), pt);
      }
      if (PdgCode == -321) {
        histReg.fill(HIST("contamination/0321in2212"), pt);
      }

      if (PdgCode == 11) {
        histReg.fill(HIST("contamination/11in2212"), pt);
      }
      if (PdgCode == -11) {
        histReg.fill(HIST("contamination/011in2212"), pt);
      }

      if (PdgCode == 13) {
        histReg.fill(HIST("contamination/13in2212"), pt);
      }
      if (PdgCode == -13) {
        histReg.fill(HIST("contamination/013in2212"), pt);
      }

    } else if (i == 2) {
      if (PdgCode == 211) {
        histReg.fill(HIST("contamination/211in321"), pt);
      }
      if (PdgCode == -211) {
        histReg.fill(HIST("contamination/0211in321"), pt);
      }

      if (PdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in321"), pt);
      }
      if (PdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in321"), pt);
      }

      if (PdgCode == 321) {
        histReg.fill(HIST("contamination/321in321"), pt);
      }
      if (PdgCode == -321) {
        histReg.fill(HIST("contamination/0321in321"), pt);
      }

      if (PdgCode == 11) {
        histReg.fill(HIST("contamination/11in321"), pt);
      }
      if (PdgCode == -11) {
        histReg.fill(HIST("contamination/011in321"), pt);
      }

      if (PdgCode == 13) {
        histReg.fill(HIST("contamination/13in321"), pt);
      }
      if (PdgCode == -13) {
        histReg.fill(HIST("contamination/013in321"), pt);
      }
    }
  }

  // nb of particles (5 particles - Pi, Pr, Ka, e, mu, and 5 antiparticles)
  static const int numParticles = 10;

  static constexpr std::string_view pidTrueRegistryNames[numParticles] = {"pidTrue/211", "pidTrue/2212", "pidTrue/321", "pidTrue/11", "pidTrue/13", "pidTrue/0211", "pidTrue/02212", "pidTrue/0321", "pidTrue/011", "pidTrue/013"};
  static constexpr std::string_view pidFalseRegistryNames[numParticles] = {"pidFalse/211", "pidFalse/2212", "pidFalse/321", "pidFalse/11", "pidFalse/13", "pidFalse/0211", "pidFalse/02212", "pidFalse/0321", "pidFalse/011", "pidFalse/013"};

  static constexpr std::string_view TPCPidTrueRegistryNames[numParticles] = {"TPCPidTrue/211", "TPCPidTrue/2212", "TPCPidTrue/321", "TPCPidTrue/11", "TPCPidTrue/13", "TPCPidTrue/0211", "TPCPidTrue/02212", "TPCPidTrue/0321", "TPCPidTrue/011", "TPCPidTrue/013"};
  static constexpr std::string_view TPCPidFalseRegistryNames[numParticles] = {"TPCPidFalse/211", "TPCPidFalse/2212", "TPCPidFalse/321", "TPCPidFalse/11", "TPCPidFalse/13", "TPCPidFalse/0211", "TPCPidFalse/02212", "TPCPidFalse/0321", "TPCPidFalse/011", "TPCPidFalse/013"};

  static constexpr std::string_view TOFPidTrueRegistryNames[numParticles] = {"TOFPidTrue/211", "TOFPidTrue/2212", "TOFPidTrue/321", "TOFPidTrue/11", "TOFPidTrue/13", "TOFPidTrue/0211", "TOFPidTrue/02212", "TOFPidTrue/0321", "TOFPidTrue/011", "TOFPidTrue/013"};
  static constexpr std::string_view TOFPidFalseRegistryNames[numParticles] = {"TOFPidFalse/211", "TOFPidFalse/2212", "TOFPidFalse/321", "TOFPidFalse/11", "TOFPidFalse/13", "TOFPidFalse/0211", "TOFPidFalse/02212", "TOFPidFalse/0321", "TOFPidFalse/011", "TOFPidFalse/013"};

  static constexpr int PdgCodes[numParticles] = {211, 2212, 321, 11, 13, -211, -2212, -321, -11, -13};
  // charge with index i corresponds to i-th particle from PdgCodes array
  static constexpr int particleCharge[numParticles] = {1, 1, 1, -1, -1, -1, -1, -1, 1, 1};
  // momentum value when to switch from pid based only on TPC (below the value) to combination of TPC and TOF (above the value)
  // i-th momentum corresponds to the i-th particle
  static constexpr float pSwitch[numParticles] = {0.5, 0.8, 0.5, 0.5, 0.5, 0.5, 0.8, 0.5, 0.5, 0.5};

  static const int maxP = 5;
  // nb of bins for TH1 hists
  static const int binsNb = 100;
  // nb of bins for TH2 hists
  static const int binsNb2D = 1000;

  HistogramRegistry histReg{
    "allHistograms",
    {{"MC/211", "MC #pi^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/11", "MC e^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/321", "MC K^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/2212", "MC p;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/13", "MC #mu^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/0211", "MC #pi^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/011", "MC e^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/0321", "MC K^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/02212", "MC #bar{p};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/013", "MC #mu^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/else", "MC else;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"pidTrue/211", "PID true #pi^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/11", "PID true e^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/321", "PID true K^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/2212", "PID true p;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/13", "PID true #mu^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/0211", "PID true #pi^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/011", "PID true e^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/0321", "PID true K^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/02212", "PID true #bar{p};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/013", "PID true #mu^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"pidFalse/211", "PID false #pi^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/11", "PID false e^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/321", "PID false K^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/2212", "PID false p;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/13", "PID false #mu^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/0211", "PID false #pi^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/011", "PID false e^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/0321", "PID false K^{-};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/02212", "PID false #bar{p};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/013", "PID false #mu^{+};p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in211", "#pi^{+} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in211", "p in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in211", "K^{+} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in211", "e^{-} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in211", "#mu^{-} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in211", "#pi^{-} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in211", "#bar{p} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in211", "K^{-} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in211", "e^{+} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in211", "#mu^{+} in #pi^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in2212", "#pi^{+} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in2212", "p in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in2212", "K^{+} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in2212", "e^{-} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in2212", "#mu^{-} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in2212", "#pi^{-} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in2212", "#bar{p} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in2212", "K^{-} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in2212", "e^{+} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in2212", "#mu^{+} in p sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in321", "#pi^{+} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in321", "p in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in321", "K^{+} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in321", "e^{-} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in321", "#mu^{-} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in321", "#pi^{-} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in321", "#bar{p} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in321", "K^{-} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in321", "e^{+} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in321", "#mu^{+} in K^{+} sample;p^{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"TPCSignalPidTrue", "PID true;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCSignalPidFalse", "PID false;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},

     {"TPCPidTrue/211", "PID true #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/11", "PID true e^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/321", "PID true K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/2212", " PID true p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/13", "PID true #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/211", "PID false #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/11", "PID false e^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/321", "PID false K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/2212", "PID false p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/13", "PID false #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},

     {"TPCPidTrue/0211", "PID true #pi^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/011", "PID true e^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/0321", "PID true K^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/02212", " PID true #bar{p};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/013", "PID true #mu^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/0211", "PID false #pi^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/011", "PID false e^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/0321", "PID false K^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/02212", "PID false #bar{p};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/013", "PID false #mu^{-};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},

     {"TOFSignalPidTrue", "PID true;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFSignalPidFalse", "PID false;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},

     {"TOFPidTrue/211", "PID true #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/2212", "PID true p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/321", "PID true K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/13", "PID true #mu^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/11", "PID true e^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/211", "PID false #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/2212", "PID false p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/321", "PID false K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/13", "PID false #mu^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/11", "PID false e^{-};p (GeV/c); TOF b#eta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},

     {"TOFPidTrue/0211", "PID true #pi^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/02212", "PID true #bar{p};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/0321", "PID true K^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/013", "PID true #mu^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/011", "PID true e^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/0211", "PID false #pi^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/02212", "PID false #bar{p};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/0321", "PID false K^{-};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/013", "PID false #mu^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/011", "PID false e^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}}}};

  template <std::size_t i, typename T>
  void pidSimple(const T& track, const int PdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arrLen)
  {
    /*
    Simsplest possibile PID, accept particle when:
    TPCSignal < X if p < Value or
    sqrt(TPCSignal^2 + TOFSignal^2) < X if p > Value
    */
    const float p = track.p();

    bool isPidFalse = false;

    if ((p < pSwitch[i]) & (track.sign() == particleCharge[i])) {
      if (abs(tpcNSigmas[i]) < nsigmacut.value) {
        if (PdgCode == PdgCodes[i]) {
          histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
          histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
          histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
          histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
          histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
        } else {
          isPidFalse = true;
        }
      }
    } else if ((p >= pSwitch[i]) & (track.sign() == particleCharge[i])) {
      if (sqrt(pow(tpcNSigmas[i], 2) + pow(tofNSigmas[i], 2)) < nsigmacut.value) {
        if (PdgCode == PdgCodes[i]) {
          histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
          histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
          histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
          histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
          histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
        } else {
          isPidFalse = true;
        }
      }
    }

    if (isPidFalse) {
      histReg.fill(HIST(pidFalseRegistryNames[i]), track.pt());
      histReg.fill(HIST(TPCPidFalseRegistryNames[i]), track.p(), track.tpcSignal());
      histReg.fill(HIST("TPCSignalPidFalse"), track.p(), track.tpcSignal());
      histReg.fill(HIST("TOFSignalPidFalse"), track.p(), track.beta());
      histReg.fill(HIST(TOFPidFalseRegistryNames[i]), track.p(), track.beta());

      double pt = track.pt();
      fillContaminationRegistry(i, PdgCode, pt);
    }
  }

  template <std::size_t i, typename T>
  void pidMinStrategy(const T& track, const int PdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arrLen)
  {
    const float p = track.p();

    // list of Nsigmas for particles
    float particleNSigma[arrLen] = {};

    bool isPidFalse = false;

    // calculate Nsigmas for every particle
    for (int j = 0; j < arrLen; ++j) {
      if (p < pSwitch[j]) {
        particleNSigma[j] = abs(tpcNSigmas[j]);
      } else if (p >= pSwitch[j]) {
        particleNSigma[j] = combindeSignal(tpcNSigmas[j], tofNSigmas[j]);
      }
    }

    if ((p < pSwitch[i]) & (track.sign() == particleCharge[i])) {
      float tmp_NSigma = abs(tpcNSigmas[i]);
      if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arrLen) == i)) {
        if (PdgCode == PdgCodes[i]) {
          histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
          histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
          histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
          histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
          histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
        } else {
          isPidFalse = true;
        }
      }
    } else if ((p >= pSwitch[i]) & (track.sign() == particleCharge[i])) {
      float tmp_NSigma = combindeSignal(tpcNSigmas[i], tofNSigmas[i]);
      if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arrLen) == i)) {
        if (PdgCode == PdgCodes[i]) {
          histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
          histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
          histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
          histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
          histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
        } else {
          isPidFalse = true;
        }
      }
    }

    if (isPidFalse) {
      histReg.fill(HIST(pidFalseRegistryNames[i]), track.pt());
      histReg.fill(HIST(TPCPidFalseRegistryNames[i]), track.p(), track.tpcSignal());
      histReg.fill(HIST("TPCSignalPidFalse"), track.p(), track.tpcSignal());
      histReg.fill(HIST("TOFSignalPidFalse"), track.p(), track.beta());
      histReg.fill(HIST(TOFPidFalseRegistryNames[i]), track.p(), track.beta());

      double pt = track.pt();
      fillContaminationRegistry(i, PdgCode, pt);
    }
  }

  template <std::size_t i, typename T>
  void pidExclusiveStrategy(const T& track, const int PdgCode, const float tpcNSigmas[], const float tofNSigmas[], int arrLen)
  {
    const float p = track.p();

    bool isPidFalse = false;

    // list of Nsigmas for particles
    float particleNSigma[arrLen] = {};

    // calculate Nsigmas for every particle
    for (int j = 0; j < arrLen; ++j) {
      if (p < pSwitch[j]) {
        particleNSigma[j] = abs(tpcNSigmas[j]);
      } else if (p >= pSwitch[j]) {
        particleNSigma[j] = combindeSignal(tpcNSigmas[j], tofNSigmas[j]);
      }
    }

    // check how many particles satisfy the condition
    int counts = 0;
    for (int j = 0; j < arrLen; ++j) {
      if (particleNSigma[j] < nsigmacut.value) {
        counts++;
      }
    }

    if (counts == 1) {
      if ((p < pSwitch[i]) & (track.sign() == particleCharge[i])) {
        float tmp_NSigma = abs(tpcNSigmas[i]);
        if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arrLen) == i)) {
          if (PdgCode == PdgCodes[i]) {
            histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
            histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
            histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
            histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
            histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
          } else {
            isPidFalse = true;
          }
        }
      } else if ((p >= pSwitch[i]) & (track.sign() == particleCharge[i])) {
        float tmp_NSigma = combindeSignal(tpcNSigmas[i], tofNSigmas[i]);
        if ((tmp_NSigma < nsigmacut.value) & (indexofSmallestElement(particleNSigma, arrLen) == i)) {
          if (PdgCode == PdgCodes[i]) {
            histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
            histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
            histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
            histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
            histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
          } else {
            isPidFalse = true;
          }
        }
      }
      if (isPidFalse) {
        histReg.fill(HIST(pidFalseRegistryNames[i]), track.pt());
        histReg.fill(HIST(TPCPidFalseRegistryNames[i]), track.p(), track.tpcSignal());
        histReg.fill(HIST("TPCSignalPidFalse"), track.p(), track.tpcSignal());
        histReg.fill(HIST("TOFSignalPidFalse"), track.p(), track.beta());
        histReg.fill(HIST(TOFPidFalseRegistryNames[i]), track.p(), track.beta());

        double pt = track.pt();
        fillContaminationRegistry(i, PdgCode, pt);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    histReg.get<TH1>(HIST("MC/211"))->Sumw2();
    histReg.get<TH1>(HIST("MC/2212"))->Sumw2();
    histReg.get<TH1>(HIST("MC/321"))->Sumw2();
    histReg.get<TH1>(HIST("MC/11"))->Sumw2();
    histReg.get<TH1>(HIST("MC/13"))->Sumw2();
    histReg.get<TH1>(HIST("MC/0211"))->Sumw2();
    histReg.get<TH1>(HIST("MC/02212"))->Sumw2();
    histReg.get<TH1>(HIST("MC/0321"))->Sumw2();
    histReg.get<TH1>(HIST("MC/011"))->Sumw2();
    histReg.get<TH1>(HIST("MC/013"))->Sumw2();

    histReg.get<TH1>(HIST("pidTrue/211"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/2212"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/321"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/11"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/13"))->Sumw2();

    histReg.get<TH1>(HIST("pidTrue/0211"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/02212"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/0321"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/011"))->Sumw2();
    histReg.get<TH1>(HIST("pidTrue/013"))->Sumw2();

    histReg.get<TH1>(HIST("pidFalse/211"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/2212"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/321"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/11"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/13"))->Sumw2();

    histReg.get<TH1>(HIST("pidFalse/0211"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/02212"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/0321"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/011"))->Sumw2();
    histReg.get<TH1>(HIST("pidFalse/013"))->Sumw2();
  }

  Configurable<float> nsigmacut{"nsigmacut", 2.5, "Value of the NSigma cut"};
  Configurable<int> strategy{"Strategy", 1, "1-Simple, 2-Minimum, 3-Exclusive"};

  Filter trackFilter = aod::track::isGlobalTrack == static_cast<uint8_t>(true);
  using pidTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TracksExtended, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTPCEl, aod::pidTPCMu, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa, aod::pidTOFEl, aod::pidTOFMu>>;
  void process(pidTracks const& tracks, aod::McParticles const& mcParticles)
  {
    for (auto& track : tracks) {
      auto particle = track.mcParticle();
      int PdgCode = particle.pdgCode();

      // pions
      if (PdgCode == 211) {
        histReg.fill(HIST("MC/211"), track.pt());
      } else if (PdgCode == -211) {
        histReg.fill(HIST("MC/0211"), track.pt());
      }
      // protons
      else if (PdgCode == 2212) {
        histReg.fill(HIST("MC/2212"), track.pt());
      } else if (PdgCode == -2212) {
        histReg.fill(HIST("MC/02212"), track.pt());
      }
      // kaons
      else if (PdgCode == 321) {
        histReg.fill(HIST("MC/321"), track.pt());
      } else if (PdgCode == -321) {
        histReg.fill(HIST("MC/0321"), track.pt());
      }
      // electrons
      else if (PdgCode == 11) {
        histReg.fill(HIST("MC/11"), track.pt());
        ;
      } else if (PdgCode == -11) {
        histReg.fill(HIST("MC/011"), track.pt());
      }
      // muons
      else if (PdgCode == 13) {
        histReg.fill(HIST("MC/13"), track.pt());
      } else if (PdgCode == -13) {
        histReg.fill(HIST("MC/013"), track.pt());
      } else {
        histReg.fill(HIST("MC/else"), track.pt());
      }

      // PID
      // indexes for particles: 0-Pi, 1-Pr, 2-Ka, 3-El, 4-Mu
      const float tpcNSigmas[numParticles] = {track.tpcNSigmaPi(), track.tpcNSigmaPr(), track.tpcNSigmaKa(), track.tpcNSigmaEl(), track.tpcNSigmaMu()};
      const float tofNSigmas[numParticles] = {track.tofNSigmaPi(), track.tofNSigmaPr(), track.tofNSigmaKa(), track.tofNSigmaEl(), track.tofNSigmaMu()};

      if (strategy.value == 1) {
        // a) strategia 1. tylko cutNSigma
        pidSimple<0>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<1>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<2>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<3>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<4>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<5>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<6>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<7>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<8>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidSimple<9>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
      } else if (strategy.value == 2) {
        // b) strategia 2. NCutSigma + warunek na minimalne NSigma
        pidMinStrategy<0>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<1>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<2>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<3>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<4>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<5>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<6>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<7>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<8>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidMinStrategy<9>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
      } else if (strategy.value == 3) {
        // c) strategia 3. exclusive NSigma
        pidExclusiveStrategy<0>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<1>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<2>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<3>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<4>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<5>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<6>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<7>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<8>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
        pidExclusiveStrategy<9>(track, PdgCode, tpcNSigmas, tofNSigmas, 5);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<pid>(cfgc),
  };
}
