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
/// \file   spectraTOF.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-12-03
/// \brief  Header for the spectraTOF task for the analysis of the spectra with the TOF and TPC detectors.
///

#ifndef PWGLF_TASKS_SPECTRATOF_H_
#define PWGLF_TASKS_SPECTRATOF_H_

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

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr PID::ID Np = 9;
static constexpr PID::ID NpCharge = Np * 2;
static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "{}^{3}He", "#alpha"};
static constexpr const char* pTCharge[NpCharge] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "d", "t", "{}^{3}He", "#alpha",
                                                   "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}", "#bar{d}", "#bar{t}", "{}^{3}#bar{He}", "#bar{#alpha}"};
static constexpr int PDGs[NpCharge] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040,
                                       -kElectron, -kMuonMinus, -kPiPlus, -kKPlus, -kProton, -1000010020, -1000010030, -1000020030, -1000020040};

static constexpr std::string_view hnsigmatpctof[NpCharge] = {"nsigmatpctof/pos/el", "nsigmatpctof/pos/mu", "nsigmatpctof/pos/pi",
                                                             "nsigmatpctof/pos/ka", "nsigmatpctof/pos/pr", "nsigmatpctof/pos/de",
                                                             "nsigmatpctof/pos/tr", "nsigmatpctof/pos/he", "nsigmatpctof/pos/al",
                                                             "nsigmatpctof/neg/el", "nsigmatpctof/neg/mu", "nsigmatpctof/neg/pi",
                                                             "nsigmatpctof/neg/ka", "nsigmatpctof/neg/pr", "nsigmatpctof/neg/de",
                                                             "nsigmatpctof/neg/tr", "nsigmatpctof/neg/he", "nsigmatpctof/neg/al"};
static constexpr std::string_view hnsigmatof[NpCharge] = {"nsigmatof/pos/el", "nsigmatof/pos/mu", "nsigmatof/pos/pi",
                                                          "nsigmatof/pos/ka", "nsigmatof/pos/pr", "nsigmatof/pos/de",
                                                          "nsigmatof/pos/tr", "nsigmatof/pos/he", "nsigmatof/pos/al",
                                                          "nsigmatof/neg/el", "nsigmatof/neg/mu", "nsigmatof/neg/pi",
                                                          "nsigmatof/neg/ka", "nsigmatof/neg/pr", "nsigmatof/neg/de",
                                                          "nsigmatof/neg/tr", "nsigmatof/neg/he", "nsigmatof/neg/al"};
static constexpr std::string_view hnsigmatpc[NpCharge] = {"nsigmatpc/pos/el", "nsigmatpc/pos/mu", "nsigmatpc/pos/pi",
                                                          "nsigmatpc/pos/ka", "nsigmatpc/pos/pr", "nsigmatpc/pos/de",
                                                          "nsigmatpc/pos/tr", "nsigmatpc/pos/he", "nsigmatpc/pos/al",
                                                          "nsigmatpc/neg/el", "nsigmatpc/neg/mu", "nsigmatpc/neg/pi",
                                                          "nsigmatpc/neg/ka", "nsigmatpc/neg/pr", "nsigmatpc/neg/de",
                                                          "nsigmatpc/neg/tr", "nsigmatpc/neg/he", "nsigmatpc/neg/al"};
static constexpr std::string_view hdeltatof[NpCharge] = {"deltatof/pos/el", "deltatof/pos/mu", "deltatof/pos/pi",
                                                         "deltatof/pos/ka", "deltatof/pos/pr", "deltatof/pos/de",
                                                         "deltatof/pos/tr", "deltatof/pos/he", "deltatof/pos/al",
                                                         "deltatof/neg/el", "deltatof/neg/mu", "deltatof/neg/pi",
                                                         "deltatof/neg/ka", "deltatof/neg/pr", "deltatof/neg/de",
                                                         "deltatof/neg/tr", "deltatof/neg/he", "deltatof/neg/al"};
static constexpr std::string_view hdeltatpc[NpCharge] = {"deltatpc/pos/el", "deltatpc/pos/mu", "deltatpc/pos/pi",
                                                         "deltatpc/pos/ka", "deltatpc/pos/pr", "deltatpc/pos/de",
                                                         "deltatpc/pos/tr", "deltatpc/pos/he", "deltatpc/pos/al",
                                                         "deltatpc/neg/el", "deltatpc/neg/mu", "deltatpc/neg/pi",
                                                         "deltatpc/neg/ka", "deltatpc/neg/pr", "deltatpc/neg/de",
                                                         "deltatpc/neg/tr", "deltatpc/neg/he", "deltatpc/neg/al"};
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

//********************************************RD**********************************************************************************************
static constexpr std::string_view hpt_numtof_str[NpCharge] = {"MC/el/pos/str/pt/numtof", "MC/mu/pos/str/pt/numtof", "MC/pi/pos/str/pt/numtof",
                                                              "MC/ka/pos/str/pt/numtof", "MC/pr/pos/str/pt/numtof", "MC/de/pos/str/pt/numtof",
                                                              "MC/tr/pos/str/pt/numtof", "MC/he/pos/str/pt/numtof", "MC/al/pos/str/pt/numtof",
                                                              "MC/el/neg/str/pt/numtof", "MC/mu/neg/str/pt/numtof", "MC/pi/neg/str/pt/numtof",
                                                              "MC/ka/neg/str/pt/numtof", "MC/pr/neg/str/pt/numtof", "MC/de/neg/str/pt/numtof",
                                                              "MC/tr/neg/str/pt/numtof", "MC/he/neg/str/pt/numtof", "MC/al/neg/str/pt/numtof"};
static constexpr std::string_view hpt_numtof_mat[NpCharge] = {"MC/el/pos/mat/pt/numtof", "MC/mu/pos/mat/pt/numtof", "MC/pi/pos/mat/pt/numtof",
                                                              "MC/ka/pos/mat/pt/numtof", "MC/pr/pos/mat/pt/numtof", "MC/de/pos/mat/pt/numtof",
                                                              "MC/tr/pos/mat/pt/numtof", "MC/he/pos/mat/pt/numtof", "MC/al/pos/mat/pt/numtof",
                                                              "MC/el/neg/mat/pt/numtof", "MC/mu/neg/mat/pt/numtof", "MC/pi/neg/mat/pt/numtof",
                                                              "MC/ka/neg/mat/pt/numtof", "MC/pr/neg/mat/pt/numtof", "MC/de/neg/mat/pt/numtof",
                                                              "MC/tr/neg/mat/pt/numtof", "MC/he/neg/mat/pt/numtof", "MC/al/neg/mat/pt/numtof"};
static constexpr std::string_view hpt_den_prm[NpCharge] = {"MC/el/pos/prm/pt/den", "MC/mu/pos/prm/pt/den", "MC/pi/pos/prm/pt/den",
                                                           "MC/ka/pos/prm/pt/den", "MC/pr/pos/prm/pt/den", "MC/de/pos/prm/pt/den",
                                                           "MC/tr/pos/prm/pt/den", "MC/he/pos/prm/pt/den", "MC/al/pos/prm/pt/den",
                                                           "MC/el/neg/prm/pt/den", "MC/mu/neg/prm/pt/den", "MC/pi/neg/prm/pt/den",
                                                           "MC/ka/neg/prm/pt/den", "MC/pr/neg/prm/pt/den", "MC/de/neg/prm/pt/den",
                                                           "MC/tr/neg/prm/pt/den", "MC/he/neg/prm/pt/den", "MC/al/neg/prm/pt/den"};
static constexpr std::string_view hpt_den_prm_recoev[NpCharge] = {"MC/el/pos/prm/pt/denrecoev", "MC/mu/pos/prm/pt/denrecoev", "MC/pi/pos/prm/pt/denrecoev",
                                                                  "MC/ka/pos/prm/pt/denrecoev", "MC/pr/pos/prm/pt/denrecoev", "MC/de/pos/prm/pt/denrecoev",
                                                                  "MC/tr/pos/prm/pt/denrecoev", "MC/he/pos/prm/pt/denrecoev", "MC/al/pos/prm/pt/denrecoev",
                                                                  "MC/el/neg/prm/pt/denrecoev", "MC/mu/neg/prm/pt/denrecoev", "MC/pi/neg/prm/pt/denrecoev",
                                                                  "MC/ka/neg/prm/pt/denrecoev", "MC/pr/neg/prm/pt/denrecoev", "MC/de/neg/prm/pt/denrecoev",
                                                                  "MC/tr/neg/prm/pt/denrecoev", "MC/he/neg/prm/pt/denrecoev", "MC/al/neg/prm/pt/denrecoev"};
static constexpr std::string_view hpt_den_prm_evsel[NpCharge] = {"MC/el/pos/prm/pt/denevsel", "MC/mu/pos/prm/pt/denevsel", "MC/pi/pos/prm/pt/denevsel",
                                                                 "MC/ka/pos/prm/pt/denevsel", "MC/pr/pos/prm/pt/denevsel", "MC/de/pos/prm/pt/denevsel",
                                                                 "MC/tr/pos/prm/pt/denevsel", "MC/he/pos/prm/pt/denevsel", "MC/al/pos/prm/pt/denevsel",
                                                                 "MC/el/neg/prm/pt/denevsel", "MC/mu/neg/prm/pt/denevsel", "MC/pi/neg/prm/pt/denevsel",
                                                                 "MC/ka/neg/prm/pt/denevsel", "MC/pr/neg/prm/pt/denevsel", "MC/de/neg/prm/pt/denevsel",
                                                                 "MC/tr/neg/prm/pt/denevsel", "MC/he/neg/prm/pt/denevsel", "MC/al/neg/prm/pt/denevsel"};
static constexpr std::string_view hpt_den_prm_goodev[NpCharge] = {"MC/el/pos/prm/pt/dengoodev", "MC/mu/pos/prm/pt/dengoodev", "MC/pi/pos/prm/pt/dengoodev",
                                                                  "MC/ka/pos/prm/pt/dengoodev", "MC/pr/pos/prm/pt/dengoodev", "MC/de/pos/prm/pt/dengoodev",
                                                                  "MC/tr/pos/prm/pt/dengoodev", "MC/he/pos/prm/pt/dengoodev", "MC/al/pos/prm/pt/dengoodev",
                                                                  "MC/el/neg/prm/pt/dengoodev", "MC/mu/neg/prm/pt/dengoodev", "MC/pi/neg/prm/pt/dengoodev",
                                                                  "MC/ka/neg/prm/pt/dengoodev", "MC/pr/neg/prm/pt/dengoodev", "MC/de/neg/prm/pt/dengoodev",
                                                                  "MC/tr/neg/prm/pt/dengoodev", "MC/he/neg/prm/pt/dengoodev", "MC/al/neg/prm/pt/dengoodev"};
static constexpr std::string_view hpt_den_prm_mcgoodev[NpCharge] = {"MC/el/pos/prm/pt/denmcgoodev", "MC/mu/pos/prm/pt/denmcgoodev", "MC/pi/pos/prm/pt/denmcgoodev",
                                                                    "MC/ka/pos/prm/pt/denmcgoodev", "MC/pr/pos/prm/pt/denmcgoodev", "MC/de/pos/prm/pt/denmcgoodev",
                                                                    "MC/tr/pos/prm/pt/denmcgoodev", "MC/he/pos/prm/pt/denmcgoodev", "MC/al/pos/prm/pt/denmcgoodev",
                                                                    "MC/el/neg/prm/pt/denmcgoodev", "MC/mu/neg/prm/pt/denmcgoodev", "MC/pi/neg/prm/pt/denmcgoodev",
                                                                    "MC/ka/neg/prm/pt/denmcgoodev", "MC/pr/neg/prm/pt/denmcgoodev", "MC/de/neg/prm/pt/denmcgoodev",
                                                                    "MC/tr/neg/prm/pt/denmcgoodev", "MC/he/neg/prm/pt/denmcgoodev", "MC/al/neg/prm/pt/denmcgoodev"};
static constexpr std::string_view hpt_den_prm_mcbadev[NpCharge] = {"MC/el/pos/prm/pt/denmcbadev", "MC/mu/pos/prm/pt/denmcbadev", "MC/pi/pos/prm/pt/denmcbadev",
                                                                   "MC/ka/pos/prm/pt/denmcbadev", "MC/pr/pos/prm/pt/denmcbadev", "MC/de/pos/prm/pt/denmcbadev",
                                                                   "MC/tr/pos/prm/pt/denmcbadev", "MC/he/pos/prm/pt/denmcbadev", "MC/al/pos/prm/pt/denmcbadev",
                                                                   "MC/el/neg/prm/pt/denmcbadev", "MC/mu/neg/prm/pt/denmcbadev", "MC/pi/neg/prm/pt/denmcbadev",
                                                                   "MC/ka/neg/prm/pt/denmcbadev", "MC/pr/neg/prm/pt/denmcbadev", "MC/de/neg/prm/pt/denmcbadev",
                                                                   "MC/tr/neg/prm/pt/denmcbadev", "MC/he/neg/prm/pt/denmcbadev", "MC/al/neg/prm/pt/denmcbadev"};
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
static constexpr std::string_view hdcaxyprmgoodevs[NpCharge] = {"dcaxyprmgoodevs/pos/el", "dcaxyprmgoodevs/pos/mu", "dcaxyprmgoodevs/pos/pi",
                                                                "dcaxyprmgoodevs/pos/ka", "dcaxyprmgoodevs/pos/pr", "dcaxyprmgoodevs/pos/de",
                                                                "dcaxyprmgoodevs/pos/tr", "dcaxyprmgoodevs/pos/he", "dcaxyprmgoodevs/pos/al",
                                                                "dcaxyprmgoodevs/neg/el", "dcaxyprmgoodevs/neg/mu", "dcaxyprmgoodevs/neg/pi",
                                                                "dcaxyprmgoodevs/neg/ka", "dcaxyprmgoodevs/neg/pr", "dcaxyprmgoodevs/neg/de",
                                                                "dcaxyprmgoodevs/neg/tr", "dcaxyprmgoodevs/neg/he", "dcaxyprmgoodevs/neg/al"};
static constexpr std::string_view hdcazprm[NpCharge] = {"dcazprm/pos/el", "dcazprm/pos/mu", "dcazprm/pos/pi",
                                                        "dcazprm/pos/ka", "dcazprm/pos/pr", "dcazprm/pos/de",
                                                        "dcazprm/pos/tr", "dcazprm/pos/he", "dcazprm/pos/al",
                                                        "dcazprm/neg/el", "dcazprm/neg/mu", "dcazprm/neg/pi",
                                                        "dcazprm/neg/ka", "dcazprm/neg/pr", "dcazprm/neg/de",
                                                        "dcazprm/neg/tr", "dcazprm/neg/he", "dcazprm/neg/al"};
static constexpr std::string_view hdcazprmgoodevs[NpCharge] = {"dcazprmgoodevs/pos/el", "dcazprmgoodevs/pos/mu", "dcazprmgoodevs/pos/pi",
                                                               "dcazprmgoodevs/pos/ka", "dcazprmgoodevs/pos/pr", "dcazprmgoodevs/pos/de",
                                                               "dcazprmgoodevs/pos/tr", "dcazprmgoodevs/pos/he", "dcazprmgoodevs/pos/al",
                                                               "dcazprmgoodevs/neg/el", "dcazprmgoodevs/neg/mu", "dcazprmgoodevs/neg/pi",
                                                               "dcazprmgoodevs/neg/ka", "dcazprmgoodevs/neg/pr", "dcazprmgoodevs/neg/de",
                                                               "dcazprmgoodevs/neg/tr", "dcazprmgoodevs/neg/he", "dcazprmgoodevs/neg/al"};
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

#endif // PWGLF_TASKS_SPECTRATOF_H_
