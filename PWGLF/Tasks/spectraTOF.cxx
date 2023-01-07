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
/// \file   spectraTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task for the analysis of the spectra with the TOF detector.
///        Depending on the configuration it can also run on tiny tables.
///

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/StaticFor.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "spectraTOF.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Spectra task
struct tofSpectra {
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 3, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", 5, "Last cluster to require in TRD for track selection. -1 does not require any TRD cluster"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, "Binning of the pT axis"};
  ConfigurableAxis binsnsigmaTPC{"binsnsigmaTPC", {200, -10, 10}, "Binning of the nsigmaTPC axis"};
  ConfigurableAxis binsnsigmaTOF{"binsnsigmaTOF", {200, -10, 10}, "Binning of the nsigmaTOF axis"};
  ConfigurableAxis binsdeltaTPC{"binsdeltaTPC", {500, -1000, 1000}, "Binning of the nsigmaTPC axis"};
  ConfigurableAxis binsdeltaTOF{"binsdeltaTOF", {500, -1000, 1000}, "Binning of the nsigmaTOF axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 100}, "Binning for multiplicity"};
  ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Binning for percentiles"};
  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 0, "Flag to use a multiplicity estimator: 0 no multiplicity, 1 MultFV0M, 2 MultFT0M, 3 MultFDDM, 4 MultTracklets, 5 MultTPC, 6 MultNTracksPV, 7 MultNTracksPVeta1, 8 CentralityFT0C"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // Standard process functions
    // Full
    if (doprocessFullEl) {
      LOG(info) << "Enabling process function processFullEl";
    }
    if (doprocessFullMu) {
      LOG(info) << "Enabling process function processFullMu";
    }
    if (doprocessFullPi) {
      LOG(info) << "Enabling process function processFullPi";
    }
    if (doprocessFullKa) {
      LOG(info) << "Enabling process function processFullKa";
    }
    if (doprocessFullPr) {
      LOG(info) << "Enabling process function processFullPr";
    }
    if (doprocessFullDe) {
      LOG(info) << "Enabling process function processFullDe";
    }
    if (doprocessFullTr) {
      LOG(info) << "Enabling process function processFullTr";
    }
    if (doprocessFullHe) {
      LOG(info) << "Enabling process function processFullHe";
    }
    if (doprocessFullAl) {
      LOG(info) << "Enabling process function processFullAl";
    }
    // LF Full
    if (doprocessLfFullEl) {
      LOG(info) << "Enabling process function processLfFullEl";
    }
    if (doprocessLfFullMu) {
      LOG(info) << "Enabling process function processLfFullMu";
    }
    if (doprocessLfFullPi) {
      LOG(info) << "Enabling process function processLfFullPi";
    }
    if (doprocessLfFullKa) {
      LOG(info) << "Enabling process function processLfFullKa";
    }
    if (doprocessLfFullPr) {
      LOG(info) << "Enabling process function processLfFullPr";
    }
    if (doprocessLfFullDe) {
      LOG(info) << "Enabling process function processLfFullDe";
    }
    if (doprocessLfFullTr) {
      LOG(info) << "Enabling process function processLfFullTr";
    }
    if (doprocessLfFullHe) {
      LOG(info) << "Enabling process function processLfFullHe";
    }
    if (doprocessLfFullAl) {
      LOG(info) << "Enabling process function processLfFullAl";
    }

    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("event/vertexz", "", HistType::kTH1F, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(3, "posZ passed");

    histos.add("track/trdSignal", "", HistType::kTH2F, {pAxis, {1000, 0, 1000, "TRD signal (a.u.)"}});
    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Eta passed");
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed");

    histos.add("Centrality/FT0M", "FT0M", HistType::kTH1F, {{binsPercentile, "Centrality FT0M"}});
    histos.add("Centrality/FT0A", "FT0A", HistType::kTH1F, {{binsPercentile, "Centrality FT0A"}});
    histos.add("Centrality/FT0C", "FT0C", HistType::kTH1F, {{binsPercentile, "Centrality FT0C"}});

    histos.add("Mult/FV0M", "MultFV0M", HistType::kTH1F, {{binsMultiplicity, "MultFV0M"}});
    histos.add("Mult/FT0M", "MultFT0M", HistType::kTH1F, {{binsMultiplicity, "MultFT0M"}});
    histos.add("Mult/FDDM", "MultFDDM", HistType::kTH1F, {{binsMultiplicity, "MultFDDM"}});

    histos.add("Mult/Tracklets", "MultTracklets", HistType::kTH1F, {{binsMultiplicity, "MultTracklets"}});
    histos.add("Mult/TPC", "MultTPC", HistType::kTH1F, {{binsMultiplicity, "MultTPC"}});
    histos.add("Mult/NTracksPV", "MultNTracksPV", HistType::kTH1F, {{binsMultiplicity, "MultNTracksPV"}});
    histos.add("Mult/NTracksPVeta1", "MultNTracksPVeta1", HistType::kTH1F, {{binsMultiplicity, "MultNTracksPVeta1"}});

    const AxisSpec dcaXyAxis{600, -3.005, 2.995, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{600, -3.005, 2.995, "DCA_{z} (cm)"};

    // 4 detectors
    histos.add("Data/pos/pt/its_tpc_trd_tof", "pos ITS-TPC-TRD-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_trd_tof", "neg ITS-TPC-TRD-TOF", kTH1F, {ptAxis});

    // 3 detectors
    histos.add("Data/pos/pt/its_tpc_trd", "pos ITS-TPC-TRD", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_trd", "neg ITS-TPC-TRD", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/its_tpc_tof", "pos ITS-TPC-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tpc_tof", "neg ITS-TPC-TOF", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/tpc_trd_tof", "pos TPC-TRD-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/tpc_trd_tof", "neg TPC-TRD-TOF", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/its_trd_tof", "pos ITS-TRD-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_trd_tof", "neg ITS-TRD-TOF", kTH1F, {ptAxis});

    // 2 detectors
    histos.add("Data/pos/pt/its_tpc", "pos ITS-TPC", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tpc", "neg ITS-TPC", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/trd_tof", "pos TRD-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/trd_tof", "neg TRD-TOF", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/tpc_tof", "pos TPC-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/tpc_tof", "neg TPC-TOF", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/its_tof", "pos ITS-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tof", "neg ITS-TOF", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/tpc_trd", "pos TPC-TRD", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/tpc_trd", "neg TPC-TRD", kTH1F, {ptAxis});

    histos.add("Data/pos/pt/its_trd", "pos ITS-TRD", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_trd", "neg ITS-TRD", kTH1F, {ptAxis});

    if (doprocessMC) {
      histos.add("MC/fake/pos", "Fake positive tracks", kTH1F, {ptAxis});
      histos.add("MC/fake/neg", "Fake negative tracks", kTH1F, {ptAxis});
    }

    for (int i = 0; i < NpCharge; i++) {

      switch (i) {
        case 0:
        case Np:
          if (doprocessFullEl == false && doprocessLfFullEl == false) {
            continue;
          }
          break;
        case 1:
        case Np + 1:
          if (doprocessFullMu == false && doprocessLfFullMu == false) {
            continue;
          }
          break;
        case 2:
        case Np + 2:
          if (doprocessFullPi == false && doprocessLfFullPi == false) {
            continue;
          }
          break;
        case 3:
        case Np + 3:
          if (doprocessFullKa == false && doprocessLfFullKa == false) {
            continue;
          }
          break;
        case 4:
        case Np + 4:
          if (doprocessFullPr == false && doprocessLfFullPr == false) {
            continue;
          }
          break;
        case 5:
        case Np + 5:
          if (doprocessFullDe == false && doprocessLfFullDe == false) {
            continue;
          }
          break;
        case 6:
        case Np + 6:
          if (doprocessFullTr == false && doprocessLfFullTr == false) {
            continue;
          }
          break;
        case 7:
        case Np + 7:
          if (doprocessFullHe == false && doprocessLfFullHe == false) {
            continue;
          }
          break;
        case 8:
        case Np + 8:
          if (doprocessFullAl == false && doprocessLfFullAl == false) {
            continue;
          }
          break;
      }

      const AxisSpec nsigmaTPCAxis{binsnsigmaTPC, Form("N_{#sigma}^{TPC}(%s)", pTCharge[i])};
      const AxisSpec nsigmaTOFAxis{binsnsigmaTOF, Form("N_{#sigma}^{TOF}(%s)", pTCharge[i])};
      const AxisSpec deltaTPCAxis{binsdeltaTPC, Form("#Delta^{TPC}(%s)", pTCharge[i])};
      const AxisSpec deltaTOFAxis{binsdeltaTOF, Form("#Delta^{TOF}(%s)", pTCharge[i])};

      histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, nsigmaTOFAxis});

      switch (multiplicityEstimator) {
        case 0: // No multiplicity
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH2F, {ptAxis, nsigmaTOFAxis});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH2F, {ptAxis, nsigmaTPCAxis});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH2F, {ptAxis, deltaTOFAxis});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH2F, {ptAxis, deltaTPCAxis});
          break;
        case 1: // MultFV0M
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultFV0M"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultFV0M"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultFV0M"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultFV0M"}});
          break;
        case 2: // MultFT0M
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultFT0M"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultFT0M"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultFT0M"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultFT0M"}});
          break;
        case 3: // MultFDDM
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultFDDM"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultFDDM"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultFDDM"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultFDDM"}});
          break;
        case 4: // MultTracklets
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultTracklets"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultTracklets"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultTracklets"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultTracklets"}});
          break;
        case 5: // MultTPC
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultTPC"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultTPC"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultTPC"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultTPC"}});
          break;
        case 6: // MultNTracksPV
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultNTracksPV"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultNTracksPV"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultNTracksPV"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultNTracksPV"}});
          break;
        case 7: // MultNTracksPVeta1
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsMultiplicity, "MultNTracksPVeta1"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsMultiplicity, "MultNTracksPVeta1"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsMultiplicity, "MultNTracksPVeta1"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsMultiplicity, "MultNTracksPVeta1"}});
          break;
        case 8: // Centrality FT0C
          histos.add(hnsigmatof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTOFAxis, {binsPercentile, "Centrality FT0C"}});
          histos.add(hnsigmatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, {binsPercentile, "Centrality FT0C"}});
          histos.add(hdeltatof[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTOFAxis, {binsPercentile, "Centrality FT0C"}});
          histos.add(hdeltatpc[i].data(), pTCharge[i], kTH3F, {ptAxis, deltaTPCAxis, {binsPercentile, "Centrality FT0C"}});
          break;
        default:
          LOG(fatal) << "Unrecognized option for multiplicity " << multiplicityEstimator;
      }
      histos.add(hdcaxy[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaXyAxis});
      histos.add(hdcaz[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaZAxis});
      histos.add(hdcaxyphi[i].data(), Form("%s -- 0.9 < #it{p}_{T} < 1.1 GeV/#it{c}", pTCharge[i]), kTH2F, {phiAxis, dcaXyAxis});

      if (doprocessMC) {
        histos.add(hpt_num_prm[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_numtof_prm[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_num_str[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_num_mat[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_den_prm[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_den_str[i].data(), pTCharge[i], kTH1F, {ptAxis});
        histos.add(hpt_den_mat[i].data(), pTCharge[i], kTH1F, {ptAxis});

        histos.add(hdcaxyprm[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaXyAxis});
        histos.add(hdcazprm[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaZAxis});
        histos.add(hdcaxystr[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaXyAxis});
        histos.add(hdcazstr[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaZAxis});
        histos.add(hdcaxymat[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaXyAxis});
        histos.add(hdcazmat[i].data(), pTCharge[i], kTH2F, {ptAxis, dcaZAxis});
      }
    }
  }

  template <bool fillFullInfo, PID::ID id, typename T, typename C>
  void fillParticleHistos(const T& track, const C& collision)
  {
    if (abs(track.rapidity(PID::getMass(id))) > cfgCutY) {
      return;
    }
    const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<id>(track);
    const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<id>(track);
    // const auto id = track.sign() > 0 ? id : id + Np;
    float multiplicity = 0.f;

    switch (multiplicityEstimator) {
      case 1: // MultFV0M
        // multiplicity = collision.multFV0M();
        // multiplicity = collision.multZeqFV0A() + collision.multZeqFV0C();
        multiplicity = collision.multZeqFV0A();
        break;
      case 2: // MultFT0M
        // multiplicity = collision.multFT0M();
        multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
        break;
      case 3: // MultFDDM
        // multiplicity = collision.multFDDM();
        multiplicity = collision.multZeqFDDA() + collision.multZeqFDDC();
        break;
      case 4: // MultTracklets
        multiplicity = collision.multTracklets();
        break;
      case 5: // MultTPC
        multiplicity = collision.multTPC();
        break;
      case 6: // MultNTracksPV
        // multiplicity = collision.multNTracksPV();
        multiplicity = collision.multZeqNTracksPV();
        break;
      case 7: // MultNTracksPVeta1
        multiplicity = collision.multNTracksPVeta1();
        break;
      case 8: // Centrality FT0C
        multiplicity = collision.centFT0C();
        break;
    }

    if (multiplicityEstimator == 0) {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatpc[id]), track.pt(), nsigmaTPC);
      } else {
        histos.fill(HIST(hnsigmatpc[id + Np]), track.pt(), nsigmaTPC);
      }
    } else {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatpc[id]), track.pt(), nsigmaTPC, multiplicity);
      } else {
        histos.fill(HIST(hnsigmatpc[id + Np]), track.pt(), nsigmaTPC, multiplicity);
      }
    }

    if constexpr (fillFullInfo) {
      const auto& deltaTPC = o2::aod::pidutils::tpcExpSignalDiff<id>(track);
      if (multiplicityEstimator == 0) {
        if (track.sign() > 0) {
          histos.fill(HIST(hdeltatpc[id]), track.pt(), deltaTPC);
        } else {
          histos.fill(HIST(hdeltatpc[id + Np]), track.pt(), deltaTPC);
        }
      } else {
        if (track.sign() > 0) {
          histos.fill(HIST(hdeltatpc[id]), track.pt(), deltaTPC, multiplicity);
        } else {
          histos.fill(HIST(hdeltatpc[id + Np]), track.pt(), deltaTPC, multiplicity);
        }
      }
    }

    if (!track.hasTOF()) {
      return;
    }

    if (track.hasTRD() && lastRequiredTrdCluster > 0) {
      int lastLayer = 0;
      for (int l = 7; l >= 0; l--) {
        if (track.trdPattern() & (1 << l)) {
          lastLayer = l;
          break;
        }
      }
      if (lastLayer < lastRequiredTrdCluster) {
        return;
      }
    }

    if (multiplicityEstimator == 0) {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatof[id]), track.pt(), nsigmaTOF);
      } else {
        histos.fill(HIST(hnsigmatof[id + Np]), track.pt(), nsigmaTOF);
      }
    } else {
      if (track.sign() > 0) {
        histos.fill(HIST(hnsigmatof[id]), track.pt(), nsigmaTOF, multiplicity);
      } else {
        histos.fill(HIST(hnsigmatof[id + Np]), track.pt(), nsigmaTOF, multiplicity);
      }
    }

    if (track.sign() > 0) {
      histos.fill(HIST(hnsigmatpctof[id]), track.pt(), nsigmaTPC, nsigmaTOF);
    } else {
      histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), nsigmaTPC, nsigmaTOF);
    }

    if constexpr (fillFullInfo) {
      const auto& deltaTOF = o2::aod::pidutils::tofExpSignalDiff<id>(track);
      if (multiplicityEstimator == 0) {
        if (track.sign() > 0) {
          histos.fill(HIST(hdeltatof[id]), track.pt(), deltaTOF);
        } else {
          histos.fill(HIST(hdeltatof[id + Np]), track.pt(), deltaTOF);
        }
      } else {
        if (track.sign() > 0) {
          histos.fill(HIST(hdeltatof[id]), track.pt(), deltaTOF, multiplicity);
        } else {
          histos.fill(HIST(hdeltatof[id + Np]), track.pt(), deltaTOF, multiplicity);
        }
      }
    }

    // Filling DCA info with the TPC+TOF PID
    if (std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.f) {
      if (track.sign() > 0) {
        histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
      } else {
        histos.fill(HIST(hdcaxy[id + Np]), track.pt(), track.dcaXY());
      }
      if (track.sign() > 0) {
        histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
      } else {
        histos.fill(HIST(hdcaz[id + Np]), track.pt(), track.dcaZ());
      }
      if (track.pt() < 1.1 && track.pt() > 0.9) {
        if (track.sign() > 0) {
          histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY());
        } else {
          histos.fill(HIST(hdcaxyphi[id + Np]), track.phi(), track.dcaXY());
        }
      }
    }
    if (!track.isGlobalTrack()) {
      return;
    }

    if constexpr (fillFullInfo) {
    }
  }

  template <bool fillHistograms = false, bool fillMultiplicity = false, typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 1);
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 2);
    }
    if (abs(collision.posZ()) > cfgCutVertex) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 3);
      histos.fill(HIST("event/vertexz"), collision.posZ());

      if constexpr (fillMultiplicity) {
        // histos.fill(HIST("Centrality/FT0M"), collision.centFT0M());
        // histos.fill(HIST("Centrality/FT0A"), collision.centFT0A());
        histos.fill(HIST("Centrality/FT0C"), collision.centFT0C());

        histos.fill(HIST("Mult/FV0M"), collision.multZeqFV0A());
        histos.fill(HIST("Mult/FT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
        histos.fill(HIST("Mult/FDDM"), collision.multZeqFDDA() + collision.multZeqFDDC());

        histos.fill(HIST("Mult/Tracklets"), collision.multTracklets());
        histos.fill(HIST("Mult/TPC"), collision.multTPC());
        histos.fill(HIST("Mult/NTracksPV"), collision.multZeqNTracksPV());
        histos.fill(HIST("Mult/NTracksPVeta1"), collision.multNTracksPVeta1());
      }
    }
    return true;
  }

  template <bool fillHistograms, typename TrackType>
  bool isTrackSelected(TrackType const& track)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 1);
    }
    if (abs(track.eta()) > cfgCutEta) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 2);
    }
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 3);
    }
    if constexpr (fillHistograms) {
      if (track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_trd_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_trd"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_trd_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_trd_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTPC()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tpc"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tpc"), track.pt());
        }
      }
      if (track.hasTRD() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/trd_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/trd_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_tof"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTOF()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_tof"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_tof"), track.pt());
        }
      }
      if (track.hasTPC() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/tpc_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/tpc_trd"), track.pt());
        }
      }
      if (track.hasITS() && track.hasTRD()) {
        if (track.sign() > 0) {
          histos.fill(HIST("Data/pos/pt/its_trd"), track.pt());
        } else {
          histos.fill(HIST("Data/neg/pt/its_trd"), track.pt());
        }
      }
    }
    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  // using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidEvTimeFlags, aod::TrackSelection, aod::TOFSignal>;

  void process(CollisionCandidate::iterator const& collision,
               TrackCandidates const& tracks)
  {
    if (!isEventSelected<true, true>(collision)) {
      return;
    }
    for (const auto& track : tracks) {
      if (!isTrackSelected<true>(track)) {
        continue;
      }
    }
  } // end of the process function

#define makeProcessFunction(processorName, inputPid, particleId, isFull, tofTable, tpcTable)   \
  void process##processorName##inputPid(CollisionCandidate::iterator const& collision,         \
                                        soa::Join<TrackCandidates,                             \
                                                  aod::pid##tofTable##inputPid,                \
                                                  aod::pid##tpcTable##inputPid> const& tracks) \
  {                                                                                            \
    if (!isEventSelected<false, false>(collision)) {                                           \
      return;                                                                                  \
    }                                                                                          \
    for (const auto& track : tracks) {                                                         \
      if (!isTrackSelected<false>(track)) {                                                    \
        continue;                                                                              \
      }                                                                                        \
      fillParticleHistos<isFull, PID::particleId>(track, collision);                           \
    }                                                                                          \
  }                                                                                            \
  PROCESS_SWITCH(tofSpectra, process##processorName##inputPid, Form("Process for the %s hypothesis from %s tables", #particleId, #processorName), false);

// Full tables
#define makeProcessFunctionFull(inputPid, particleId) makeProcessFunction(Full, inputPid, particleId, true, TOFFull, TPCFull)

  makeProcessFunctionFull(El, Electron);
  makeProcessFunctionFull(Mu, Muon);
  makeProcessFunctionFull(Pi, Pion);
  makeProcessFunctionFull(Ka, Kaon);
  makeProcessFunctionFull(Pr, Proton);
  makeProcessFunctionFull(De, Deuteron);
  makeProcessFunctionFull(Tr, Triton);
  makeProcessFunctionFull(He, Helium3);
  makeProcessFunctionFull(Al, Alpha);
#undef makeProcessFunctionFull

// Full LF tables
#define makeProcessFunctionFull(inputPid, particleId) makeProcessFunction(LfFull, inputPid, particleId, true, TOFFull, TPCLfFull)

  makeProcessFunctionFull(El, Electron);
  makeProcessFunctionFull(Mu, Muon);
  makeProcessFunctionFull(Pi, Pion);
  makeProcessFunctionFull(Ka, Kaon);
  makeProcessFunctionFull(Pr, Proton);
  makeProcessFunctionFull(De, Deuteron);
  makeProcessFunctionFull(Tr, Triton);
  makeProcessFunctionFull(He, Helium3);
  makeProcessFunctionFull(Al, Alpha);
#undef makeProcessFunctionFull

  template <std::size_t i, typename T1, typename T2>
  void fillHistograms_MC(T1 const& tracks, T2 const& mcParticles)
  {

    switch (i) {
      case 0:
      case Np:
        if (doprocessFullEl == false && doprocessLfFullEl == false) {
          return;
        }
        break;
      case 1:
      case Np + 1:
        if (doprocessFullMu == false && doprocessLfFullMu == false) {
          return;
        }
        break;
      case 2:
      case Np + 2:
        if (doprocessFullPi == false && doprocessLfFullPi == false) {
          return;
        }
        break;
      case 3:
      case Np + 3:
        if (doprocessFullKa == false && doprocessLfFullKa == false) {
          return;
        }
        break;
      case 4:
      case Np + 4:
        if (doprocessFullPr == false && doprocessLfFullPr == false) {
          return;
        }
        break;
      case 5:
      case Np + 5:
        if (doprocessFullDe == false && doprocessLfFullDe == false) {
          return;
        }
        break;
      case 6:
      case Np + 6:
        if (doprocessFullTr == false && doprocessLfFullTr == false) {
          return;
        }
        break;
      case 7:
      case Np + 7:
        if (doprocessFullHe == false && doprocessLfFullHe == false) {
          return;
        }
        break;
      case 8:
      case Np + 8:
        if (doprocessFullAl == false && doprocessLfFullAl == false) {
          return;
        }
        break;
    }

    for (auto& track : tracks) {
      if (!track.isGlobalTrackWoDCA()) {
        continue;
      }
      if (!track.has_mcParticle()) {
        if (track.sign() > 0) {
          histos.fill(HIST("MC/fake/pos"), track.pt());
        } else {
          histos.fill(HIST("MC/fake/neg"), track.pt());
        }
        continue;
      }
      const auto mcParticle = track.mcParticle();

      if (mcParticle.pdgCode() != PDGs[i]) {
        continue;
      }
      if (std::abs(mcParticle.eta()) > cfgCutEta) {
        continue;
      }
      if (std::abs(mcParticle.y()) > cfgCutY) {
        continue;
      }
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST(hdcaxystr[i]), track.pt(), track.dcaXY());
          histos.fill(HIST(hdcazstr[i]), track.pt(), track.dcaZ());
        } else {
          histos.fill(HIST(hdcaxymat[i]), track.pt(), track.dcaXY());
          histos.fill(HIST(hdcazmat[i]), track.pt(), track.dcaZ());
        }
      } else {
        histos.fill(HIST(hdcaxyprm[i]), track.pt(), track.dcaXY());
        histos.fill(HIST(hdcazprm[i]), track.pt(), track.dcaZ());
      }

      if (!track.isGlobalTrack()) { // Skipping tracks that don't pass the standard cuts
        continue;
      }

      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == 4) {
          histos.fill(HIST(hpt_num_str[i]), track.pt());
        } else {
          histos.fill(HIST(hpt_num_mat[i]), track.pt());
        }
      } else {
        histos.fill(HIST(hpt_num_prm[i]), track.pt());

        if (track.hasTRD() && lastRequiredTrdCluster > 0) {
          int lastLayer = 0;
          for (int l = 7; l >= 0; l--) {
            if (track.trdPattern() & (1 << l)) {
              lastLayer = l;
              break;
            }
          }
          if (lastLayer < lastRequiredTrdCluster) {
            continue;
          }
        }

        if (track.hasTOF()) {
          histos.fill(HIST(hpt_numtof_prm[i]), track.pt());
        }
      }
    }

    for (auto& particle : mcParticles) {
      if (std::abs(particle.eta()) > cfgCutEta) {
        continue;
      }
      if (std::abs(particle.y()) > cfgCutY) {
        continue;
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (particle.pdgCode() != PDGs[i]) {
        continue;
      }

      if (!particle.isPhysicalPrimary()) {
        if (particle.getProcess() == 4) {
          histos.fill(HIST(hpt_den_str[i]), particle.pt());
        } else {
          histos.fill(HIST(hpt_den_mat[i]), particle.pt());
        }
      } else {
        histos.fill(HIST(hpt_den_prm[i]), particle.pt());
      }
    }
  }

  void processMC(soa::Join<aod::Tracks, aod::TracksExtra,
                           aod::TracksDCA, aod::McTrackLabels,
                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::TrackSelection> const& tracks,
                 const aod::McParticles& mcParticles)
  {
    // LOGF(info, "Enter processMC!");
    static_for<0, 17>([&](auto i) {
      fillHistograms_MC<i>(tracks, mcParticles);
    });
  }
  PROCESS_SWITCH(tofSpectra, processMC, "Process MC", false);

}; // end of spectra task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofSpectra>(cfgc)};
}
