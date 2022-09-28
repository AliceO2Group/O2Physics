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
#include "Common/Core/TrackSelection.h"
#include "Framework/StaticFor.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Spectra task
struct tofSpectra {
  static constexpr PID::ID Np = 9;
  static constexpr PID::ID NpCharge = Np * 2;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "{}^{3}He", "#alpha"};
  static constexpr const char* pTCharge[NpCharge] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "d", "t", "{}^{3}He", "#alpha",
                                                     "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}", "#bar{d}", "#bar{t}", "{}^{3}#bar{He}", "#bar{#alpha}"};
  static constexpr int PDGs[NpCharge] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040,
                                         -kElectron, -kMuonMinus, -kPiPlus, -kKPlus, -kProton, -1000010020, -1000010030, -1000020030, -1000020040};
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 3, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0}, ""};
  Configurable<bool> isRun2{"isRun2", false, "Flag to process Run 2 data"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr std::string_view hp[NpCharge] = {"p/pos/el", "p/pos/mu", "p/pos/pi",
                                                    "p/pos/ka", "p/pos/pr", "p/pos/de",
                                                    "p/pos/tr", "p/pos/he", "p/pos/al",
                                                    "p/neg/el", "p/neg/mu", "p/neg/pi",
                                                    "p/neg/ka", "p/neg/pr", "p/neg/de",
                                                    "p/neg/tr", "p/neg/he", "p/neg/al"};
  static constexpr std::string_view hpt[NpCharge] = {"pt/pos/el", "pt/pos/mu", "pt/pos/pi",
                                                     "pt/pos/ka", "pt/pos/pr", "pt/pos/de",
                                                     "pt/pos/tr", "pt/pos/he", "pt/pos/al",
                                                     "pt/neg/el", "pt/neg/mu", "pt/neg/pi",
                                                     "pt/neg/ka", "pt/neg/pr", "pt/neg/de",
                                                     "pt/neg/tr", "pt/neg/he", "pt/neg/al"};
  static constexpr std::string_view hnsigmatpctof[NpCharge] = {"nsigmatpctof/pos/el", "nsigmatpctof/pos/mu", "nsigmatpctof/pos/pi",
                                                               "nsigmatpctof/pos/ka", "nsigmatpctof/pos/pr", "nsigmatpctof/pos/de",
                                                               "nsigmatpctof/pos/tr", "nsigmatpctof/pos/he", "nsigmatpctof/pos/al",
                                                               "nsigmatpctof/neg/el", "nsigmatpctof/neg/mu", "nsigmatpctof/neg/pi",
                                                               "nsigmatpctof/neg/ka", "nsigmatpctof/neg/pr", "nsigmatpctof/neg/de",
                                                               "nsigmatpctof/neg/tr", "nsigmatpctof/neg/he", "nsigmatpctof/neg/al"};
  static constexpr std::string_view hdcaxy[NpCharge] = {"dcaxy/pos/el", "dcaxy/pos/mu", "dcaxy/pos/pi",
                                                        "dcaxy/pos/ka", "dcaxy/pos/pr", "dcaxy/pos/de",
                                                        "dcaxy/pos/tr", "dcaxy/pos/he", "dcaxy/pos/al",
                                                        "dcaxy/neg/el", "dcaxy/neg/mu", "dcaxy/neg/pi",
                                                        "dcaxy/neg/ka", "dcaxy/neg/pr", "dcaxy/neg/de",
                                                        "dcaxy/neg/tr", "dcaxy/neg/he", "dcaxy/neg/al"};
  static constexpr std::string_view hdcaz[NpCharge] = {"dcaz/pos/el", "dcaz/pos/mu", "dcaz/pos/pi",
                                                       "dcaz/pos/ka", "dcaz/pos/pr", "dcaz/pos/de",
                                                       "dcaz/pos/tr", "dcaz/pos/he", "dcaz/pos/al",
                                                       "dcaz/neg/el", "dcaz/neg/mu", "dcaz/neg/pi",
                                                       "dcaz/neg/ka", "dcaz/neg/pr", "dcaz/neg/de",
                                                       "dcaz/neg/tr", "dcaz/neg/he", "dcaz/neg/al"};
  static constexpr std::string_view hdcaxyphi[NpCharge] = {"dcaxyphi/pos/el", "dcaxyphi/pos/mu", "dcaxyphi/pos/pi",
                                                           "dcaxyphi/pos/ka", "dcaxyphi/pos/pr", "dcaxyphi/pos/de",
                                                           "dcaxyphi/pos/tr", "dcaxyphi/pos/he", "dcaxyphi/pos/al",
                                                           "dcaxyphi/neg/el", "dcaxyphi/neg/mu", "dcaxyphi/neg/pi",
                                                           "dcaxyphi/neg/ka", "dcaxyphi/neg/pr", "dcaxyphi/neg/de",
                                                           "dcaxyphi/neg/tr", "dcaxyphi/neg/he", "dcaxyphi/neg/al"};
  // MC
  static constexpr std::string_view hpt_num_prm[NpCharge] = {"MC/el/pos/prm/pt/num", "MC/mu/pos/prm/pt/num", "MC/pi/pos/prm/pt/num",
                                                             "MC/ka/pos/prm/pt/num", "MC/pr/pos/prm/pt/num", "MC/de/pos/prm/pt/num",
                                                             "MC/tr/pos/prm/pt/num", "MC/he/pos/prm/pt/num", "MC/al/pos/prm/pt/num",
                                                             "MC/el/neg/prm/pt/num", "MC/mu/neg/prm/pt/num", "MC/pi/neg/prm/pt/num",
                                                             "MC/ka/neg/prm/pt/num", "MC/pr/neg/prm/pt/num", "MC/de/neg/prm/pt/num",
                                                             "MC/tr/neg/prm/pt/num", "MC/he/neg/prm/pt/num", "MC/al/neg/prm/pt/num"};
  static constexpr std::string_view hpt_numtof_prm[NpCharge] = {"MC/el/pos/prm/pt/numtof", "MC/mu/pos/prm/pt/numtof", "MC/pi/pos/prm/pt/numtof",
                                                                "MC/ka/pos/prm/pt/numtof", "MC/pr/pos/prm/pt/numtof", "MC/de/pos/prm/pt/numtof",
                                                                "MC/tr/pos/prm/pt/numtof", "MC/he/pos/prm/pt/numtof", "MC/al/pos/prm/pt/numtof",
                                                                "MC/el/neg/prm/pt/numtof", "MC/mu/neg/prm/pt/numtof", "MC/pi/neg/prm/pt/numtof",
                                                                "MC/ka/neg/prm/pt/numtof", "MC/pr/neg/prm/pt/numtof", "MC/de/neg/prm/pt/numtof",
                                                                "MC/tr/neg/prm/pt/numtof", "MC/he/neg/prm/pt/numtof", "MC/al/neg/prm/pt/numtof"};
  static constexpr std::string_view hpt_den_prm[NpCharge] = {"MC/el/pos/prm/pt/den", "MC/mu/pos/prm/pt/den", "MC/pi/pos/prm/pt/den",
                                                             "MC/ka/pos/prm/pt/den", "MC/pr/pos/prm/pt/den", "MC/de/pos/prm/pt/den",
                                                             "MC/tr/pos/prm/pt/den", "MC/he/pos/prm/pt/den", "MC/al/pos/prm/pt/den",
                                                             "MC/el/neg/prm/pt/den", "MC/mu/neg/prm/pt/den", "MC/pi/neg/prm/pt/den",
                                                             "MC/ka/neg/prm/pt/den", "MC/pr/neg/prm/pt/den", "MC/de/neg/prm/pt/den",
                                                             "MC/tr/neg/prm/pt/den", "MC/he/neg/prm/pt/den", "MC/al/neg/prm/pt/den"};
  static constexpr std::string_view hpt_num_str[NpCharge] = {"MC/el/pos/str/pt/num", "MC/mu/pos/str/pt/num", "MC/pi/pos/str/pt/num",
                                                             "MC/ka/pos/str/pt/num", "MC/pr/pos/str/pt/num", "MC/de/pos/str/pt/num",
                                                             "MC/tr/pos/str/pt/num", "MC/he/pos/str/pt/num", "MC/al/pos/str/pt/num",
                                                             "MC/el/neg/str/pt/num", "MC/mu/neg/str/pt/num", "MC/pi/neg/str/pt/num",
                                                             "MC/ka/neg/str/pt/num", "MC/pr/neg/str/pt/num", "MC/de/neg/str/pt/num",
                                                             "MC/tr/neg/str/pt/num", "MC/he/neg/str/pt/num", "MC/al/neg/str/pt/num"};
  static constexpr std::string_view hpt_den_str[NpCharge] = {"MC/el/pos/str/pt/den", "MC/mu/pos/str/pt/den", "MC/pi/pos/str/pt/den",
                                                             "MC/ka/pos/str/pt/den", "MC/pr/pos/str/pt/den", "MC/de/pos/str/pt/den",
                                                             "MC/tr/pos/str/pt/den", "MC/he/pos/str/pt/den", "MC/al/pos/str/pt/den",
                                                             "MC/el/neg/str/pt/den", "MC/mu/neg/str/pt/den", "MC/pi/neg/str/pt/den",
                                                             "MC/ka/neg/str/pt/den", "MC/pr/neg/str/pt/den", "MC/de/neg/str/pt/den",
                                                             "MC/tr/neg/str/pt/den", "MC/he/neg/str/pt/den", "MC/al/neg/str/pt/den"};
  static constexpr std::string_view hpt_num_mat[NpCharge] = {"MC/el/pos/mat/pt/num", "MC/mu/pos/mat/pt/num", "MC/pi/pos/mat/pt/num",
                                                             "MC/ka/pos/mat/pt/num", "MC/pr/pos/mat/pt/num", "MC/de/pos/mat/pt/num",
                                                             "MC/tr/pos/mat/pt/num", "MC/he/pos/mat/pt/num", "MC/al/pos/mat/pt/num",
                                                             "MC/el/neg/mat/pt/num", "MC/mu/neg/mat/pt/num", "MC/pi/neg/mat/pt/num",
                                                             "MC/ka/neg/mat/pt/num", "MC/pr/neg/mat/pt/num", "MC/de/neg/mat/pt/num",
                                                             "MC/tr/neg/mat/pt/num", "MC/he/neg/mat/pt/num", "MC/al/neg/mat/pt/num"};
  static constexpr std::string_view hpt_den_mat[NpCharge] = {"MC/el/pos/mat/pt/den", "MC/mu/pos/mat/pt/den", "MC/pi/pos/mat/pt/den",
                                                             "MC/ka/pos/mat/pt/den", "MC/pr/pos/mat/pt/den", "MC/de/pos/mat/pt/den",
                                                             "MC/tr/pos/mat/pt/den", "MC/he/pos/mat/pt/den", "MC/al/pos/mat/pt/den",
                                                             "MC/el/neg/mat/pt/den", "MC/mu/neg/mat/pt/den", "MC/pi/neg/mat/pt/den",
                                                             "MC/ka/neg/mat/pt/den", "MC/pr/neg/mat/pt/den", "MC/de/neg/mat/pt/den",
                                                             "MC/tr/neg/mat/pt/den", "MC/he/neg/mat/pt/den", "MC/al/neg/mat/pt/den"};
  static constexpr std::string_view hdcaxyprm[NpCharge] = {"dcaxyprm/pos/el", "dcaxyprm/pos/mu", "dcaxyprm/pos/pi",
                                                           "dcaxyprm/pos/ka", "dcaxyprm/pos/pr", "dcaxyprm/pos/de",
                                                           "dcaxyprm/pos/tr", "dcaxyprm/pos/he", "dcaxyprm/pos/al",
                                                           "dcaxyprm/neg/el", "dcaxyprm/neg/mu", "dcaxyprm/neg/pi",
                                                           "dcaxyprm/neg/ka", "dcaxyprm/neg/pr", "dcaxyprm/neg/de",
                                                           "dcaxyprm/neg/tr", "dcaxyprm/neg/he", "dcaxyprm/neg/al"};
  static constexpr std::string_view hdcazprm[NpCharge] = {"dcazprm/pos/el", "dcazprm/pos/mu", "dcazprm/pos/pi",
                                                          "dcazprm/pos/ka", "dcazprm/pos/pr", "dcazprm/pos/de",
                                                          "dcazprm/pos/tr", "dcazprm/pos/he", "dcazprm/pos/al",
                                                          "dcazprm/neg/el", "dcazprm/neg/mu", "dcazprm/neg/pi",
                                                          "dcazprm/neg/ka", "dcazprm/neg/pr", "dcazprm/neg/de",
                                                          "dcazprm/neg/tr", "dcazprm/neg/he", "dcazprm/neg/al"};
  static constexpr std::string_view hdcaxystr[NpCharge] = {"dcaxystr/pos/el", "dcaxystr/pos/mu", "dcaxystr/pos/pi",
                                                           "dcaxystr/pos/ka", "dcaxystr/pos/pr", "dcaxystr/pos/de",
                                                           "dcaxystr/pos/tr", "dcaxystr/pos/he", "dcaxystr/pos/al",
                                                           "dcaxystr/neg/el", "dcaxystr/neg/mu", "dcaxystr/neg/pi",
                                                           "dcaxystr/neg/ka", "dcaxystr/neg/pr", "dcaxystr/neg/de",
                                                           "dcaxystr/neg/tr", "dcaxystr/neg/he", "dcaxystr/neg/al"};
  static constexpr std::string_view hdcazstr[NpCharge] = {"dcazstr/pos/el", "dcazstr/pos/mu", "dcazstr/pos/pi",
                                                          "dcazstr/pos/ka", "dcazstr/pos/pr", "dcazstr/pos/de",
                                                          "dcazstr/pos/tr", "dcazstr/pos/he", "dcazstr/pos/al",
                                                          "dcazstr/neg/el", "dcazstr/neg/mu", "dcazstr/neg/pi",
                                                          "dcazstr/neg/ka", "dcazstr/neg/pr", "dcazstr/neg/de",
                                                          "dcazstr/neg/tr", "dcazstr/neg/he", "dcazstr/neg/al"};
  static constexpr std::string_view hdcaxymat[NpCharge] = {"dcaxymat/pos/el", "dcaxymat/pos/mu", "dcaxymat/pos/pi",
                                                           "dcaxymat/pos/ka", "dcaxymat/pos/pr", "dcaxymat/pos/de",
                                                           "dcaxymat/pos/tr", "dcaxymat/pos/he", "dcaxymat/pos/al",
                                                           "dcaxymat/neg/el", "dcaxymat/neg/mu", "dcaxymat/neg/pi",
                                                           "dcaxymat/neg/ka", "dcaxymat/neg/pr", "dcaxymat/neg/de",
                                                           "dcaxymat/neg/tr", "dcaxymat/neg/he", "dcaxymat/neg/al"};
  static constexpr std::string_view hdcazmat[NpCharge] = {"dcazmat/pos/el", "dcazmat/pos/mu", "dcazmat/pos/pi",
                                                          "dcazmat/pos/ka", "dcazmat/pos/pr", "dcazmat/pos/de",
                                                          "dcazmat/pos/tr", "dcazmat/pos/he", "dcazmat/pos/al",
                                                          "dcazmat/neg/el", "dcazmat/neg/mu", "dcazmat/neg/pi",
                                                          "dcazmat/neg/ka", "dcazmat/neg/pr", "dcazmat/neg/de",
                                                          "dcazmat/neg/tr", "dcazmat/neg/he", "dcazmat/neg/al"};

  void init(o2::framework::InitContext&)
  {
    if (doprocessFullEl == true && doprocessTinyEl == true) {
      LOGF(fatal, "Cannot enable processFullEl and processTinyEl at the same time. Please choose one.");
    }
    if (doprocessFullMu == true && doprocessTinyMu == true) {
      LOGF(fatal, "Cannot enable processFullMu and processTinyMu at the same time. Please choose one.");
    }
    if (doprocessFullPi == true && doprocessTinyPi == true) {
      LOGF(fatal, "Cannot enable processFullPi and processTinyPi at the same time. Please choose one.");
    }
    if (doprocessFullKa == true && doprocessTinyKa == true) {
      LOGF(fatal, "Cannot enable processFullKa and processTinyKa at the same time. Please choose one.");
    }
    if (doprocessFullPr == true && doprocessTinyPr == true) {
      LOGF(fatal, "Cannot enable processFullPr and processTinyPr at the same time. Please choose one.");
    }
    if (doprocessFullDe == true && doprocessTinyDe == true) {
      LOGF(fatal, "Cannot enable processFullDe and processTinyDe at the same time. Please choose one.");
    }
    if (doprocessFullTr == true && doprocessTinyTr == true) {
      LOGF(fatal, "Cannot enable processFullTr and processTinyTr at the same time. Please choose one.");
    }
    if (doprocessFullHe == true && doprocessTinyHe == true) {
      LOGF(fatal, "Cannot enable processFullHe and processTinyHe at the same time. Please choose one.");
    }
    if (doprocessFullAl == true && doprocessTinyAl == true) {
      LOGF(fatal, "Cannot enable processFullAl and processTinyAl at the same time. Please choose one.");
    }

    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("event/vertexz", "", HistType::kTH1F, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(3, "posZ passed");
    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Eta passed");
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed");
    histos.add("p/Unselected", "Unselected", kTH1F, {pAxis});
    histos.add("pt/Unselected", "Unselected", kTH1F, {ptAxis});
    for (int i = 0; i < NpCharge; i++) {
      histos.add(hp[i].data(), pTCharge[i], kTH1F, {pAxis});
      histos.add(hpt[i].data(), pTCharge[i], kTH1F, {ptAxis});
    }
    // histos.add("electronbeta/hp_El", "", kTH1F, {pAxis});
    // histos.add("electronbeta/hpt_El", "", kTH1F, {ptAxis});
    // histos.add("electronbeta/hlength_El", ";Track Length (cm);Tracks", kTH1D, {{100, 0, 1000}});
    // histos.add("electronbeta/htime_El", ";TOF Time (ns);Tracks", kTH1D, {{1000, 0, 600}});
    // histos.add("electronbeta/hp_beta_El", ";#it{p} (GeV/#it{c});#beta - #beta_{e};Tracks", kTH2D, {pAxis, {100, -0.01, 0.01}});
    // histos.add("electronbeta/hp_betasigma_El", ";#it{p} (GeV/#it{c});(#beta - #beta_{e})/#sigma;Tracks", kTH2D, {pAxis, {100, -5, 5}});

    const AxisSpec dcaXyAxis{600, -3.005, 2.995, "DCA_{xy} (cm)"};
    const AxisSpec phiAxis{200, 0, 7, "#it{#varphi} (rad)"};
    const AxisSpec dcaZAxis{600, -3.005, 2.995, "DCA_{z} (cm)"};

    histos.add("Data/pos/pt/its_tpc_tof", "pos ITS-TPC-TOF", kTH1F, {ptAxis});
    histos.add("Data/pos/pt/its_tpc", "pos ITS-TPC", kTH1F, {ptAxis});
    histos.add("Data/pos/pt/tpc", "pos TPC", kTH1F, {ptAxis});
    histos.add("Data/pos/pt/its", "pos ITS", kTH1F, {ptAxis});

    histos.add("Data/neg/pt/its_tpc_tof", "neg ITS-TPC-TOF", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its_tpc", "neg ITS-TPC", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/tpc", "neg TPC", kTH1F, {ptAxis});
    histos.add("Data/neg/pt/its", "neg ITS", kTH1F, {ptAxis});

    for (int i = 0; i < NpCharge; i++) {

      switch (i) {
        case 0:
        case Np:
          if (doprocessFullEl == false && doprocessTinyEl == false) {
            continue;
          }
          break;
        case 1:
        case Np + 1:
          if (doprocessFullMu == false && doprocessTinyMu == false) {
            continue;
          }
          break;
        case 2:
        case Np + 2:
          if (doprocessFullPi == false && doprocessTinyPi == false) {
            continue;
          }
          break;
        case 3:
        case Np + 3:
          if (doprocessFullKa == false && doprocessTinyKa == false) {
            continue;
          }
          break;
        case 4:
        case Np + 4:
          if (doprocessFullPr == false && doprocessTinyPr == false) {
            continue;
          }
          break;
        case 5:
        case Np + 5:
          if (doprocessFullDe == false && doprocessTinyDe == false) {
            continue;
          }
          break;
        case 6:
        case Np + 6:
          if (doprocessFullTr == false && doprocessTinyTr == false) {
            continue;
          }
          break;
        case 7:
        case Np + 7:
          if (doprocessFullHe == false && doprocessTinyHe == false) {
            continue;
          }
          break;
        case 8:
        case Np + 8:
          if (doprocessFullAl == false && doprocessTinyAl == false) {
            continue;
          }
          break;
      }

      const AxisSpec nsigmaTPCAxis{200, -10, 10, Form("N_{#sigma}^{TPC}(%s)", pTCharge[i])};
      const AxisSpec nsigmaTOFAxis{200, -10, 10, Form("N_{#sigma}^{TOF}(%s)", pTCharge[i])};
      histos.add(hnsigmatpctof[i].data(), pTCharge[i], kTH3F, {ptAxis, nsigmaTPCAxis, nsigmaTOFAxis});
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

  template <bool fillFullInfo, PID::ID id, typename T>
  void fillParticleHistos(const T& track)
  {
    if (abs(track.rapidity(PID::getMass(id))) > cfgCutY) {
      return;
    }
    const auto& nsigmaTOF = o2::aod::pidutils::tofNSigma<id>(track);
    const auto& nsigmaTPC = o2::aod::pidutils::tpcNSigma<id>(track);
    if (track.sign() > 0) {
      histos.fill(HIST(hnsigmatpctof[id]), track.pt(), nsigmaTPC, nsigmaTOF);
    } else {
      histos.fill(HIST(hnsigmatpctof[id + Np]), track.pt(), nsigmaTPC, nsigmaTOF);
    }

    // if (std::abs(nsigmaTOF) < 2) {
    if (std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.f) {
      if (track.sign() > 0) {
        histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
        histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
        if (track.pt() < 1.1 && track.pt() > 0.9) {
          histos.fill(HIST(hdcaxyphi[id]), track.phi(), track.dcaXY());
        }
      } else {
        histos.fill(HIST(hdcaxy[id + Np]), track.pt(), track.dcaXY());
        histos.fill(HIST(hdcaz[id + Np]), track.pt(), track.dcaZ());
        if (track.pt() < 1.1 && track.pt() > 0.9) {
          histos.fill(HIST(hdcaxyphi[id + Np]), track.phi(), track.dcaXY());
        }
      }
    }
    if (!track.isGlobalTrack()) {
      return;
    }
    if (abs(nsigmaTOF) > cfgNSigmaCut) {
      return;
    }
    if (track.sign() > 0) {
      histos.fill(HIST(hp[id]), track.p());
      histos.fill(HIST(hpt[id]), track.pt());
    } else {
      histos.fill(HIST(hp[id + Np]), track.p());
      histos.fill(HIST(hpt[id + Np]), track.pt());
    }

    if constexpr (fillFullInfo) {
    }
  }

  template <bool fillHistograms, typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if constexpr (fillHistograms) {
      histos.fill(HIST("evsel"), 1);
    }
    if (isRun2 && !collision.sel7()) {
      return false;

    } else if (!collision.sel8()) {
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
      if (track.sign() > 0) {
        if (track.hasITS() && track.hasTPC()) {
          if (track.hasTOF()) {
            histos.fill(HIST("Data/pos/pt/its_tpc_tof"), track.pt());
          } else {
            histos.fill(HIST("Data/pos/pt/its_tpc"), track.pt());
          }
        }
        if (track.hasTPC()) {
          histos.fill(HIST("Data/pos/pt/tpc"), track.pt());
        }
        if (track.hasITS()) {
          histos.fill(HIST("Data/pos/pt/its"), track.pt());
        }
      } else {
        if (track.hasITS() && track.hasTPC()) {
          if (track.hasTOF()) {
            histos.fill(HIST("Data/neg/pt/its_tpc_tof"), track.pt());
          } else {
            histos.fill(HIST("Data/neg/pt/its_tpc"), track.pt());
          }
        }
        if (track.hasTPC()) {
          histos.fill(HIST("Data/neg/pt/tpc"), track.pt());
        }
        if (track.hasITS()) {
          histos.fill(HIST("Data/neg/pt/its"), track.pt());
        }
      }
    }
    if (!track.hasTOF()) {
      return false;
    }
    if constexpr (fillHistograms) {
      histos.fill(HIST("tracksel"), 4);
      histos.fill(HIST("p/Unselected"), track.p());
      histos.fill(HIST("pt/Unselected"), track.pt());
    }

    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                    aod::pidEvTimeFlags, aod::TrackSelection, aod::TOFSignal>;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               TrackCandidates const& tracks)
  {
    if (!isEventSelected<true>(collision)) {
      return;
    }
    for (const auto& track : tracks) {
      if (!isTrackSelected<true>(track)) {
        continue;
      }

      //
      // if (TMath::Abs(track.separationbetael() < 1.f)) {
      //   histos.fill(HIST("electronbeta/hp_El"), track.p());
      //   histos.fill(HIST("electronbeta/hpt_El"), track.pt());
      //   histos.fill(HIST("electronbeta/hlength_El"), track.length());
      //   histos.fill(HIST("electronbeta/htime_El"), track.tofSignal() / 1000);
      //   histos.fill(HIST("electronbeta/hp_beta_El"), track.p(), track.diffbetael());
      //   histos.fill(HIST("electronbeta/hp_betasigma_El"), track.p(), track.separationbetael());
      // }
    }
  } // end of the process function

// Full tables
#define makeProcessFunction(inputPid, particleId)                                \
  void processFull##inputPid(CollisionCandidate const& collision,                \
                             soa::Join<TrackCandidates,                          \
                                       aod::pidTOFFull##inputPid,                \
                                       aod::pidTPCFull##inputPid> const& tracks) \
  {                                                                              \
    if (!isEventSelected<false>(collision)) {                                    \
      return;                                                                    \
    }                                                                            \
    for (const auto& track : tracks) {                                           \
      if (!isTrackSelected<false>(track)) {                                      \
        continue;                                                                \
      }                                                                          \
      fillParticleHistos<true, PID::particleId>(track);                          \
    }                                                                            \
  }                                                                              \
  PROCESS_SWITCH(tofSpectra, processFull##inputPid, Form("Process for the %s hypothesis from full tables", #particleId), false);

  makeProcessFunction(El, Electron);
  makeProcessFunction(Mu, Muon);
  makeProcessFunction(Pi, Pion);
  makeProcessFunction(Ka, Kaon);
  makeProcessFunction(Pr, Proton);
  makeProcessFunction(De, Deuteron);
  makeProcessFunction(Tr, Triton);
  makeProcessFunction(He, Helium3);
  makeProcessFunction(Al, Alpha);
#undef makeProcessFunction

// Tiny tables
#define makeProcessFunction(inputPid, particleId)                            \
  void processTiny##inputPid(CollisionCandidate const& collision,            \
                             soa::Join<TrackCandidates,                      \
                                       aod::pidTOF##inputPid,                \
                                       aod::pidTPC##inputPid> const& tracks) \
  {                                                                          \
    if (!isEventSelected<false>(collision)) {                                \
      return;                                                                \
    }                                                                        \
    for (const auto& track : tracks) {                                       \
      if (!isTrackSelected<false>(track)) {                                  \
        continue;                                                            \
      }                                                                      \
      fillParticleHistos<true, PID::particleId>(track);                      \
    }                                                                        \
  }                                                                          \
  PROCESS_SWITCH(tofSpectra, processTiny##inputPid, Form("Process for the %s hypothesis from tiny tables", #particleId), false);

  makeProcessFunction(El, Electron);
  makeProcessFunction(Mu, Muon);
  makeProcessFunction(Pi, Pion);
  makeProcessFunction(Ka, Kaon);
  makeProcessFunction(Pr, Proton);
  makeProcessFunction(De, Deuteron);
  makeProcessFunction(Tr, Triton);
  makeProcessFunction(He, Helium3);
  makeProcessFunction(Al, Alpha);
#undef makeProcessFunction

  template <std::size_t i, typename T1, typename T2>
  void fillHistograms_MC(T1 const& tracks, T2 const& mcParticles)
  {
    for (auto& track : tracks) {
      if (!track.has_mcParticle()) {
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
      if (!track.isGlobalTrackWoDCA()) {
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
