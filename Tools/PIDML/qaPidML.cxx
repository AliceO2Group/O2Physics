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
/// \brief Task to check ML PID efficiency. Based on Maja's simpleApplyPidOnnxModel.cxx code
/// \author Łukasz Sawicki
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include <TParameter.h>
#include "Tools/PIDML/pidOnnxModel.h"
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct pidml {
  static const int maxP = 5;
  // nb of bins for TH1 hists
  static const int binsNb = 100;
  // nb of bins for TH2 hists
  static const int binsNb2D = 1000;
  static const int numParticles = 3;

  static constexpr std::string_view pidTrueRegistryNames[numParticles] = {"pidTrue/211", "pidTrue/2212", "pidTrue/321"};
  static constexpr std::string_view pidFalseRegistryNames[numParticles] = {"pidFalse/211", "pidFalse/2212", "pidFalse/321"};

  static constexpr std::string_view TPCPidTrueRegistryNames[numParticles] = {"TPCPidTrue/211", "TPCPidTrue/2212", "TPCPidTrue/321"};
  static constexpr std::string_view TPCPidFalseRegistryNames[numParticles] = {"TPCPidFalse/211", "TPCPidFalse/2212", "TPCPidFalse/321"};

  static constexpr std::string_view TOFPidTrueRegistryNames[numParticles] = {"TOFPidTrue/211", "TOFPidTrue/2212", "TOFPidTrue/321"};
  static constexpr std::string_view TOFPidFalseRegistryNames[numParticles] = {"TOFPidFalse/211", "TOFPidFalse/2212", "TOFPidFalse/321"};

  // available particles: 211, 2212, 321
  static constexpr int particlesPdgCode[numParticles] = {211, 2212, 321};

  // values of track momentum when to switch from only TPC signal to combined TPC and TOF signal
  // i-th momentum corresponds to the i-th particle
  static constexpr float pSwitchValue[numParticles] = {0.5, 0.8, 0.5};

  HistogramRegistry histReg{
    "allHistograms",
    {{"MC/211", "MC #pi^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/11", "MC e^{-};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/321", "MC K^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/2212", "MC p;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/13", "MC #mu^{-};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/0211", "MC #pi^{-};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/011", "MC e^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/0321", "MC K^{-};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/02212", "MC #bar{p};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/013", "MC #mu^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"MC/else", "MC else;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"pidTrue/211", "PID true #pi^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/321", "PID true K^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidTrue/2212", "PID true p;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"pidFalse/211", "PID false #pi^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/321", "PID false K^{+};p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"pidFalse/2212", "PID false p;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in211", "#pi^{+} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in211", "p in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in211", "K^{+} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in211", "e^{-} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in211", "#mu^{-} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in211", "#pi^{-} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in211", "#bar{p} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in211", "K^{-} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in211", "e^{+} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in211", "#mu^{+} in #pi^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in2212", "#pi^{+} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in2212", "p in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in2212", "K^{+} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in2212", "e^{-} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in2212", "#mu^{-} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in2212", "#pi^{-} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in2212", "#bar{p} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in2212", "K^{-} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in2212", "e^{+} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in2212", "#mu^{+} in p sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/211in321", "#pi^{+} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/2212in321", "p in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/321in321", "K^{+} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/11in321", "e^{-} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/13in321", "#mu^{-} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"contamination/0211in321", "#pi^{-} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/02212in321", "#bar{p} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/0321in321", "K^{-} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/011in321", "e^{+} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},
     {"contamination/013in321", "#mu^{+} in K^{+} sample;p_{T} (GeV/c);Counts", {HistType::kTH1F, {{binsNb, 0, maxP}}}},

     {"TPCSignalPidTrue", "PID true;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCSignalPidFalse", "PID false;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},

     {"TPCPidTrue/211", "PID true #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/321", "PID true K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidTrue/2212", " PID true p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/211", "PID false #pi^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/321", "PID false K^{+};p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},
     {"TPCPidFalse/2212", "PID false p;p (GeV/c); TPC Signal (dE/dx)", {HistType::kTH2F, {{binsNb2D, 0, 10}, {binsNb2D, 0, 300}}}},

     {"TOFSignalPidTrue", "PID true;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFSignalPidFalse", "PID false;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},

     {"TOFPidTrue/211", "PID true #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/2212", "PID true p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidTrue/321", "PID true K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/211", "PID false #pi^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/2212", "PID false p;p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}},
     {"TOFPidFalse/321", "PID false K^{+};p (GeV/c); TOF #beta", {HistType::kTH2F, {{binsNb2D, 0.2, 10}, {110, 0, 1.1}}}}}};

  void fillContaminationRegistry(int i, int pdgCode, double pt)
  {
    if (i == 0) {
      if (pdgCode == 211) {
        histReg.fill(HIST("contamination/211in211"), pt);
      }
      if (pdgCode == -211) {
        histReg.fill(HIST("contamination/0211in211"), pt);
      }

      if (pdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in211"), pt);
      }
      if (pdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in211"), pt);
      }

      if (pdgCode == 321) {
        histReg.fill(HIST("contamination/321in211"), pt);
      }
      if (pdgCode == -321) {
        histReg.fill(HIST("contamination/0321in211"), pt);
      }

      if (pdgCode == 11) {
        histReg.fill(HIST("contamination/11in211"), pt);
      }
      if (pdgCode == -11) {
        histReg.fill(HIST("contamination/011in211"), pt);
      }

      if (pdgCode == 13) {
        histReg.fill(HIST("contamination/13in211"), pt);
      }
      if (pdgCode == -13) {
        histReg.fill(HIST("contamination/013in211"), pt);
      }

    } else if (i == 1) {
      if (pdgCode == 211) {
        histReg.fill(HIST("contamination/211in2212"), pt);
      }
      if (pdgCode == -211) {
        histReg.fill(HIST("contamination/0211in2212"), pt);
      }

      if (pdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in2212"), pt);
      }
      if (pdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in2212"), pt);
      }

      if (pdgCode == 321) {
        histReg.fill(HIST("contamination/321in2212"), pt);
      }
      if (pdgCode == -321) {
        histReg.fill(HIST("contamination/0321in2212"), pt);
      }

      if (pdgCode == 11) {
        histReg.fill(HIST("contamination/11in2212"), pt);
      }
      if (pdgCode == -11) {
        histReg.fill(HIST("contamination/011in2212"), pt);
      }

      if (pdgCode == 13) {
        histReg.fill(HIST("contamination/13in2212"), pt);
      }
      if (pdgCode == -13) {
        histReg.fill(HIST("contamination/013in2212"), pt);
      }

    } else if (i == 2) {
      if (pdgCode == 211) {
        histReg.fill(HIST("contamination/211in321"), pt);
      }
      if (pdgCode == -211) {
        histReg.fill(HIST("contamination/0211in321"), pt);
      }

      if (pdgCode == 2212) {
        histReg.fill(HIST("contamination/2212in321"), pt);
      }
      if (pdgCode == -2212) {
        histReg.fill(HIST("contamination/02212in321"), pt);
      }

      if (pdgCode == 321) {
        histReg.fill(HIST("contamination/321in321"), pt);
      }
      if (pdgCode == -321) {
        histReg.fill(HIST("contamination/0321in321"), pt);
      }

      if (pdgCode == 11) {
        histReg.fill(HIST("contamination/11in321"), pt);
      }
      if (pdgCode == -11) {
        histReg.fill(HIST("contamination/011in321"), pt);
      }

      if (pdgCode == 13) {
        histReg.fill(HIST("contamination/13in321"), pt);
      }
      if (pdgCode == -13) {
        histReg.fill(HIST("contamination/013in321"), pt);
      }
    }
  }

  template <typename T>
  void fillMcHistos(const T& track, const int pdgCode)
  {
    // pions
    if (pdgCode == 211) {
      histReg.fill(HIST("MC/211"), track.pt());
    } else if (pdgCode == -211) {
      histReg.fill(HIST("MC/0211"), track.pt());
    }
    // protons
    else if (pdgCode == 2212) {
      histReg.fill(HIST("MC/2212"), track.pt());
    } else if (pdgCode == -2212) {
      histReg.fill(HIST("MC/02212"), track.pt());
    }
    // kaons
    else if (pdgCode == 321) {
      histReg.fill(HIST("MC/321"), track.pt());
    } else if (pdgCode == -321) {
      histReg.fill(HIST("MC/0321"), track.pt());
    }
    // electrons
    else if (pdgCode == 11) {
      histReg.fill(HIST("MC/11"), track.pt());
    } else if (pdgCode == -11) {
      histReg.fill(HIST("MC/011"), track.pt());
    }
    // muons
    else if (pdgCode == 13) {
      histReg.fill(HIST("MC/13"), track.pt());
    } else if (pdgCode == -13) {
      histReg.fill(HIST("MC/013"), track.pt());
    } else {
      histReg.fill(HIST("MC/else"), track.pt());
    }
  }

  template <std::size_t i, typename T>
  void fillPidHistos(const T& track, const int pdgCode, bool isPidTrue)
  {
    if (isPidTrue) {
      histReg.fill(HIST(pidTrueRegistryNames[i]), track.pt());
      histReg.fill(HIST(TPCPidTrueRegistryNames[i]), track.p(), track.tpcSignal());
      histReg.fill(HIST("TPCSignalPidTrue"), track.p(), track.tpcSignal());
      histReg.fill(HIST("TOFSignalPidTrue"), track.p(), track.beta());
      histReg.fill(HIST(TOFPidTrueRegistryNames[i]), track.p(), track.beta());
    } else {
      histReg.fill(HIST(pidFalseRegistryNames[i]), track.pt());
      histReg.fill(HIST(TPCPidFalseRegistryNames[i]), track.p(), track.tpcSignal());
      histReg.fill(HIST("TPCSignalPidFalse"), track.p(), track.tpcSignal());
      histReg.fill(HIST("TOFSignalPidFalse"), track.p(), track.beta());
      histReg.fill(HIST(TOFPidFalseRegistryNames[i]), track.p(), track.beta());

      double pt = track.pt();
      fillContaminationRegistry(i, pdgCode, pt);
    }
  }

  int getParticlePdg(float pidCertainties[])
  {
    // index of the biggest value in an array
    int index = 0;
    // index of the second biggest value in an array
    int smaller_index = 0;

    for (int j = 1; j < numParticles; j++) {
      if (pidCertainties[j] > pidCertainties[index]) {
        // assign new indexes
        smaller_index = index;
        index = j;
      }
    }

    // return 0 if certainty with index 'index' is below 0.5 or two indexes have the same value, else map index to particle pdgCode
    if ((pidCertainties[index] < 0.5f) | ((pidCertainties[index] == pidCertainties[smaller_index]) & (smaller_index != 0))) {
      return 0;
    } else {
      return particlesPdgCode[index];
    }
  }

  template <std::size_t i, typename T>
  void pidML(const T& track, const int pdgCodeMC)
  {
    float pidCertainties[3];
    if (track.p() < pSwitchValue[i]) {
      pidCertainties[0] = model211TPC.applyModel(track);
      pidCertainties[1] = model2212TPC.applyModel(track);
      pidCertainties[2] = model321TPC.applyModel(track);
    } else {
      pidCertainties[0] = model211All.applyModel(track);
      pidCertainties[1] = model2212All.applyModel(track);
      pidCertainties[2] = model321All.applyModel(track);
    }
    int pid = getParticlePdg(pidCertainties);
    // condition for sign: we want to work only with pi, p and K, without antiparticles
    if (pid == particlesPdgCode[i] && track.sign() == 1) {
      if (pdgCodeMC == particlesPdgCode[i]) {
        fillPidHistos<i>(track, pdgCodeMC, true);
      } else {
        fillPidHistos<i>(track, pdgCodeMC, false);
      }
    }
  }

  // one model for one particle; Model with all TPC and TOF signal
  PidONNXModel model211All;
  PidONNXModel model2212All;
  PidONNXModel model321All;
  // Model with only TPC signal model
  PidONNXModel model211TPC;
  PidONNXModel model2212TPC;
  PidONNXModel model321TPC;

  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"useCCDB", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};

  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  void init(InitContext const&)
  {
    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    } else {
      model211All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 211, kTPCTOF, 0.5f);
      model2212All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 2211, kTPCTOF, 0.5f);
      model321All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 321, kTPCTOF, 0.5f);

      model211TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 211, kTPCOnly, 0.5f);
      model2212TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 2211, kTPCOnly, 0.5f);
      model321TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, 321, kTPCOnly, 0.5f);
    }
  }

  Filter trackFilter = requireGlobalTrackInFilter();
  using pidTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::TOFSignal>>;
  void process(aod::Collisions const& collisions, pidTracks const& tracks, aod::McParticles const& /*mcParticles*/, aod::BCsWithTimestamps const&)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      model211All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 211, kTPCTOF, 0.5f);
      model2212All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 2211, kTPCTOF, 0.5f);
      model321All = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 321, kTPCTOF, 0.5f);

      model211TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 211, kTPCOnly, 0.5f);
      model2212TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 2211, kTPCOnly, 0.5f);
      model321TPC = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, bc.timestamp(), 321, kTPCOnly, 0.5f);
    }

    for (auto& track : tracks) {
      auto particle = track.mcParticle_as<aod::McParticles_000>();
      int pdgCodeMC = particle.pdgCode();

      fillMcHistos(track, pdgCodeMC);

      // only 3 particles can be predicted by model
      static_for<0, 2>([&](auto i) {
        pidML<i>(track, pdgCodeMC);
      });
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<pidml>(cfgc),
  };
}
