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

#ifndef PWGLF_DATAMODEL_SPECTRATOF_H_
#define PWGLF_DATAMODEL_SPECTRATOF_H_

#include <memory>

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"

#include "TPDGCode.h"

static constexpr o2::track::PID::ID Np = 9;
static constexpr int NCharges = 2;
static constexpr o2::track::PID::ID NpCharge = Np * NCharges;
static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "{}^{3}He", "#alpha"};
static constexpr const char* pN[Np] = {"el", "mu", "pi", "ka", "pr", "de", "tr", "he", "al"};
static constexpr const char* cN[NCharges] = {"pos", "neg"};
static constexpr const char* pTCharge[NpCharge] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "d", "t", "{}^{3}He", "#alpha",
                                                   "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}", "#bar{d}", "#bar{t}", "{}^{3}#bar{He}", "#bar{#alpha}"};
static constexpr const char* pNCharge[NpCharge] = {"pos/el", "pos/mu", "pos/pi", "pos/ka", "pos/pr", "pos/de", "pos/tr", "pos/he", "pos/al",
                                                   "neg/el", "neg/mu", "neg/pi", "neg/ka", "neg/pr", "neg/de", "neg/tr", "neg/he", "neg/al"};
static constexpr int PDGs[NpCharge] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040,
                                       -kElectron, -kMuonMinus, -kPiPlus, -kKPlus, -kProton, -1000010020, -1000010030, -1000020030, -1000020040};

std::shared_ptr<TH2> hMultiplicityvsPercentile;
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
static constexpr std::string_view hpt_mism_its_prm[NpCharge] = {"MC/el/pos/prm/pt/mismITS", "MC/mu/pos/prm/pt/mismITS", "MC/pi/pos/prm/pt/mismITS",
                                                                "MC/ka/pos/prm/pt/mismITS", "MC/pr/pos/prm/pt/mismITS", "MC/de/pos/prm/pt/mismITS",
                                                                "MC/tr/pos/prm/pt/mismITS", "MC/he/pos/prm/pt/mismITS", "MC/al/pos/prm/pt/mismITS",
                                                                "MC/el/neg/prm/pt/mismITS", "MC/mu/neg/prm/pt/mismITS", "MC/pi/neg/prm/pt/mismITS",
                                                                "MC/ka/neg/prm/pt/mismITS", "MC/pr/neg/prm/pt/mismITS", "MC/de/neg/prm/pt/mismITS",
                                                                "MC/tr/neg/prm/pt/mismITS", "MC/he/neg/prm/pt/mismITS", "MC/al/neg/prm/pt/mismITS"};
static constexpr std::string_view hpt_mism_tpc_prm[NpCharge] = {"MC/el/pos/prm/pt/mismTPC", "MC/mu/pos/prm/pt/mismTPC", "MC/pi/pos/prm/pt/mismTPC",
                                                                "MC/ka/pos/prm/pt/mismTPC", "MC/pr/pos/prm/pt/mismTPC", "MC/de/pos/prm/pt/mismTPC",
                                                                "MC/tr/pos/prm/pt/mismTPC", "MC/he/pos/prm/pt/mismTPC", "MC/al/pos/prm/pt/mismTPC",
                                                                "MC/el/neg/prm/pt/mismTPC", "MC/mu/neg/prm/pt/mismTPC", "MC/pi/neg/prm/pt/mismTPC",
                                                                "MC/ka/neg/prm/pt/mismTPC", "MC/pr/neg/prm/pt/mismTPC", "MC/de/neg/prm/pt/mismTPC",
                                                                "MC/tr/neg/prm/pt/mismTPC", "MC/he/neg/prm/pt/mismTPC", "MC/al/neg/prm/pt/mismTPC"};
static constexpr std::string_view hpt_mism_trd_prm[NpCharge] = {"MC/el/pos/prm/pt/mismTRD", "MC/mu/pos/prm/pt/mismTRD", "MC/pi/pos/prm/pt/mismTRD",
                                                                "MC/ka/pos/prm/pt/mismTRD", "MC/pr/pos/prm/pt/mismTRD", "MC/de/pos/prm/pt/mismTRD",
                                                                "MC/tr/pos/prm/pt/mismTRD", "MC/he/pos/prm/pt/mismTRD", "MC/al/pos/prm/pt/mismTRD",
                                                                "MC/el/neg/prm/pt/mismTRD", "MC/mu/neg/prm/pt/mismTRD", "MC/pi/neg/prm/pt/mismTRD",
                                                                "MC/ka/neg/prm/pt/mismTRD", "MC/pr/neg/prm/pt/mismTRD", "MC/de/neg/prm/pt/mismTRD",
                                                                "MC/tr/neg/prm/pt/mismTRD", "MC/he/neg/prm/pt/mismTRD", "MC/al/neg/prm/pt/mismTRD"};
static constexpr std::string_view hpt_mism_tof_prm[NpCharge] = {"MC/el/pos/prm/pt/mismTOF", "MC/mu/pos/prm/pt/mismTOF", "MC/pi/pos/prm/pt/mismTOF",
                                                                "MC/ka/pos/prm/pt/mismTOF", "MC/pr/pos/prm/pt/mismTOF", "MC/de/pos/prm/pt/mismTOF",
                                                                "MC/tr/pos/prm/pt/mismTOF", "MC/he/pos/prm/pt/mismTOF", "MC/al/pos/prm/pt/mismTOF",
                                                                "MC/el/neg/prm/pt/mismTOF", "MC/mu/neg/prm/pt/mismTOF", "MC/pi/neg/prm/pt/mismTOF",
                                                                "MC/ka/neg/prm/pt/mismTOF", "MC/pr/neg/prm/pt/mismTOF", "MC/de/neg/prm/pt/mismTOF",
                                                                "MC/tr/neg/prm/pt/mismTOF", "MC/he/neg/prm/pt/mismTOF", "MC/al/neg/prm/pt/mismTOF"};
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
static constexpr std::string_view hpt_numtofgoodmatch_prm[NpCharge] = {"MC/el/pos/prm/pt/numtofgoodmatch", "MC/mu/pos/prm/pt/numtofgoodmatch", "MC/pi/pos/prm/pt/numtofgoodmatch",
                                                                       "MC/ka/pos/prm/pt/numtofgoodmatch", "MC/pr/pos/prm/pt/numtofgoodmatch", "MC/de/pos/prm/pt/numtofgoodmatch",
                                                                       "MC/tr/pos/prm/pt/numtofgoodmatch", "MC/he/pos/prm/pt/numtofgoodmatch", "MC/al/pos/prm/pt/numtofgoodmatch",
                                                                       "MC/el/neg/prm/pt/numtofgoodmatch", "MC/mu/neg/prm/pt/numtofgoodmatch", "MC/pi/neg/prm/pt/numtofgoodmatch",
                                                                       "MC/ka/neg/prm/pt/numtofgoodmatch", "MC/pr/neg/prm/pt/numtofgoodmatch", "MC/de/neg/prm/pt/numtofgoodmatch",
                                                                       "MC/tr/neg/prm/pt/numtofgoodmatch", "MC/he/neg/prm/pt/numtofgoodmatch", "MC/al/neg/prm/pt/numtofgoodmatch"};

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

// Derived data model for cut variation
namespace o2::aod
{
namespace spectra
{

template <typename binningType>
typename binningType::binned_t packInTable(const float& valueToBin)
{
  if (valueToBin <= binningType::binned_min) {
    return (binningType::underflowBin);
  } else if (valueToBin >= binningType::binned_max) {
    return (binningType::overflowBin);
  } else if (valueToBin >= 0) {
    return (static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) + 0.5f));
  } else {
    return (static_cast<typename binningType::binned_t>((valueToBin / binningType::bin_width) - 0.5f));
  }
}
// Function to unpack a binned value into a float
template <typename binningType>
float unPack(const typename binningType::binned_t& valueToUnpack)
{
  return binningType::bin_width * static_cast<float>(valueToUnpack);
}

struct binningDCA {
 public:
  typedef int16_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 6.0;
  static constexpr float binned_min = -6.0;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

struct binningNSigma {
 public:
  typedef int16_t binned_t;
  static constexpr int nbins = (1 << 8 * sizeof(binned_t)) - 2;
  static constexpr binned_t overflowBin = nbins >> 1;
  static constexpr binned_t underflowBin = -(nbins >> 1);
  static constexpr float binned_max = 10.0;
  static constexpr float binned_min = -10.0;
  static constexpr float bin_width = (binned_max - binned_min) / nbins;
};

// Collision info
DECLARE_SOA_INDEX_COLUMN(BC, bc); //! Most probably BC to where this collision has occurred
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(Sel8, sel8, bool);
DECLARE_SOA_COLUMN(MultNTracksPVeta1, multNTracksPVeta1, int);
DECLARE_SOA_DYNAMIC_COLUMN(CentFV0A, centFV0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0A, centFT0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(CentFT0C, centFT0C, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFV0A, multZeqFV0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFT0A, multZeqFT0A, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFT0C, multZeqFT0C, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFDDA, multZeqFDDA, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqFDDC, multZeqFDDC, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultTracklets, multTracklets, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(MultTPC, multTPC, //! Dummy
                           [](bool /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(SelectionBit, selection_bit, //! Dummy
                           [](aod::evsel::EventSelectionFlags /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt0, isInelGt0, //! is INEL > 0
                           [](int multPveta1) -> bool { return multPveta1 > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsInelGt1, isInelGt1, //! is INEL > 1
                           [](int multPveta1) -> bool { return multPveta1 > 1; });

// Track info
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                  //! Index to the collision
DECLARE_SOA_COLUMN(PtSigned, ptSigned, float);                                   //! Pt (signed) of the track
DECLARE_SOA_COLUMN(Eta, eta, float);                                             //! Eta of the track
DECLARE_SOA_COLUMN(Phi, phi, float);                                             //! Phi of the track
DECLARE_SOA_COLUMN(EvTimeT0AC, evTimeT0AC, float);                               //! Event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(EvTimeT0ACErr, evTimeT0ACErr, float);                         //! Resolution of the event time of the track computed with the T0AC
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);                      //! IsPVContributor
DECLARE_SOA_COLUMN(DetectorMap, detectorMap, uint8_t);                           //! Detector map: see enum DetectorMapEnum
DECLARE_SOA_COLUMN(LastTRDCluster, lastTRDCluster, int8_t);                      //! Index of the last cluster in the TRD, -1 if no TRD information
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                                        //! Has or not the TRD match
DECLARE_SOA_COLUMN(TPCNSigmaStorePi, tpcNSigmaStorePi, binningNSigma::binned_t); //! Stored binned nsigma with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCNSigmaStoreKa, tpcNSigmaStoreKa, binningNSigma::binned_t); //! Stored binned nsigma with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCNSigmaStorePr, tpcNSigmaStorePr, binningNSigma::binned_t); //! Stored binned nsigma with the TPC detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaStorePi, tofNSigmaStorePi, binningNSigma::binned_t); //! Stored binned nsigma with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaStoreKa, tofNSigmaStoreKa, binningNSigma::binned_t); //! Stored binned nsigma with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaStorePr, tofNSigmaStorePr, binningNSigma::binned_t); //! Stored binned nsigma with the TOF detector for proton
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPi, tpcNSigmaPi,                             //! Unpacked NSigma TPC Pi
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaKa, tpcNSigmaKa, //! Unpacked NSigma TPC Ka
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNSigmaPr, tpcNSigmaPr, //! Unpacked NSigma TPC Pr
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPi, tofNSigmaPi, //! Unpacked NSigma TOF Pi
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaKa, tofNSigmaKa, //! Unpacked NSigma TOF Ka
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(TOFNSigmaPr, tofNSigmaPr, //! Unpacked NSigma TOF Pr
                           [](binningNSigma::binned_t binned) -> float { return unPack<binningNSigma>(binned); });

DECLARE_SOA_COLUMN(DCAxyStore, dcaxyStore, binningDCA::binned_t); //! Stored binned dcaxy
DECLARE_SOA_COLUMN(DCAzStore, dcazStore, binningDCA::binned_t);   //! Stored binned dcaz
DECLARE_SOA_DYNAMIC_COLUMN(DCAxy, dcaXY,                          //! Unpacked dcaxy
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAz, dcaZ, //! Unpacked dcaz
                           [](binningDCA::binned_t binned) -> float { return unPack<binningDCA>(binned); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Absolute value of signed pT
                           [](float signedPt) -> float { return std::abs(signedPt); });
DECLARE_SOA_DYNAMIC_COLUMN(HasITS, hasITS, //! Dummy
                           [](float /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Dummy
                           [](float /*v*/) -> bool { return true; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](float tofSignal) -> bool { return tofSignal > 0; });
DECLARE_SOA_DYNAMIC_COLUMN(TRDSignal, trdSignal, //! Dummy
                           [](float /*v*/) -> float { return 0.f; });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float signedpt, float eta) -> float { return std::abs(signedpt) * cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackType, trackType, [](float /*v*/) -> uint8_t { return o2::aod::track::TrackTypeEnum::Track; });
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                                       // if a track passed the isGlobalTrack requirement
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);                             // if a track passed the isGlobalTrackWoDCA requirement
DECLARE_SOA_DYNAMIC_COLUMN(Flags, flags, [](float /*v*/) -> uint32_t { return 0; });          // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(TRDPattern, trdPattern, [](float /*v*/) -> uint8_t { return 0; }); // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity,                                                //! Track rapidity, computed under the mass assumption given as input
                           [](float signedPt, float eta, float mass) -> float {
                             const auto pt = std::abs(signedPt);
                             const auto p = std::abs(signedPt) * cosh(eta);
                             const auto pz = std::sqrt(p * p - pt * pt);
                             const auto energy = sqrt(p * p + mass * mass);
                             return 0.5f * log((energy + pz) / (energy - pz));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(IsInAcceptanceTrack, isInAcceptanceTrack, [](float /*v*/) -> bool { return false; }); // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackITS, isQualityTrackITS, [](float /*v*/) -> bool { return false; });     // Dummy
DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrackTPC, isQualityTrackTPC, [](float /*v*/) -> bool { return false; });     // Dummy

} // namespace spectra

DECLARE_SOA_TABLE(SpColls, "AOD", "SPCOLLS",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  spectra::CentFT0M,
                  spectra::Sel8,
                  spectra::MultNTracksPVeta1,
                  spectra::RunNumber,
                  spectra::IsInelGt0<spectra::MultNTracksPVeta1>,
                  spectra::IsInelGt1<spectra::MultNTracksPVeta1>,
                  spectra::CentFV0A<spectra::Sel8>,
                  spectra::CentFT0A<spectra::Sel8>,
                  spectra::CentFT0C<spectra::Sel8>,
                  spectra::MultZeqFV0A<spectra::Sel8>,
                  spectra::MultZeqFT0A<spectra::Sel8>,
                  spectra::MultZeqFT0C<spectra::Sel8>,
                  spectra::MultZeqFDDA<spectra::Sel8>,
                  spectra::MultZeqFDDC<spectra::Sel8>,
                  spectra::MultZeqNTracksPV<spectra::Sel8>,
                  spectra::MultTracklets<spectra::Sel8>,
                  spectra::MultTPC<spectra::Sel8>,
                  spectra::SelectionBit<>);
using SpColl = SpColls::iterator;

DECLARE_SOA_TABLE(SpTracks, "AOD", "SPTRACKS",
                  o2::soa::Index<>,
                  spectra::CollisionId,
                  spectra::TPCNSigmaStorePi, spectra::TPCNSigmaStoreKa, spectra::TPCNSigmaStorePr,
                  spectra::TOFNSigmaStorePi, spectra::TOFNSigmaStoreKa, spectra::TOFNSigmaStorePr,
                  spectra::PtSigned, spectra::Eta, spectra::Phi,
                  track::Length,
                  track::TPCSignal,
                  track::TPCChi2NCl, track::ITSChi2NCl, track::TOFChi2,
                  track::TPCNClsShared,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  spectra::IsPVContributor,
                  track::ITSClusterSizes,
                  spectra::HasTRD,
                  //   pidtofevtime::EvTimeTOF,
                  //   pidtofevtime::EvTimeTOFErr,
                  //   pidtofevtime::EvTimeTOFMult,
                  // spectra::EvTimeT0AC,
                  // spectra::EvTimeT0ACErr,
                  // collision::CollisionTime,
                  // collision::CollisionTimeRes,
                  pidflags::TOFFlags,
                  spectra::DCAxyStore,
                  spectra::DCAzStore,
                  spectra::IsGlobalTrack,
                  spectra::IsGlobalTrackWoDCA,
                  spectra::DCAxy<spectra::DCAxyStore>,
                  spectra::DCAz<spectra::DCAzStore>,
                  spectra::Pt<spectra::PtSigned>,
                  track::Sign<spectra::PtSigned>,
                  spectra::P<spectra::PtSigned, spectra::Eta>,
                  spectra::Rapidity<spectra::PtSigned, spectra::Eta>,
                  spectra::HasITS<track::ITSClusterSizes>,
                  spectra::HasTPC<track::TPCChi2NCl>,
                  spectra::HasTOF<track::TOFChi2>,
                  spectra::TRDSignal<track::TOFChi2>,
                  spectra::Flags<track::TOFChi2>,
                  spectra::TrackType<track::TOFChi2>,
                  spectra::TRDPattern<track::TOFChi2>,
                  spectra::IsInAcceptanceTrack<track::TOFChi2>, // Dummy
                  spectra::IsQualityTrackITS<track::TOFChi2>,   // Dummy
                  spectra::IsQualityTrackTPC<track::TOFChi2>,   // Dummy
                  track::v001::ITSClusterMap<track::ITSClusterSizes>, track::v001::ITSNCls<track::ITSClusterSizes>, track::v001::ITSNClsInnerBarrel<track::ITSClusterSizes>,
                  track::TPCFractionSharedCls<track::TPCNClsShared, track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  pidflags::IsEvTimeDefined<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOF<pidflags::TOFFlags>,
                  pidflags::IsEvTimeT0AC<pidflags::TOFFlags>,
                  pidflags::IsEvTimeTOFT0AC<pidflags::TOFFlags>,
                  spectra::TOFNSigmaPi<spectra::TOFNSigmaStorePi>,
                  spectra::TOFNSigmaKa<spectra::TOFNSigmaStoreKa>,
                  spectra::TOFNSigmaPr<spectra::TOFNSigmaStorePr>,
                  spectra::TPCNSigmaPi<spectra::TPCNSigmaStorePi>,
                  spectra::TPCNSigmaKa<spectra::TPCNSigmaStoreKa>,
                  spectra::TPCNSigmaPr<spectra::TPCNSigmaStorePr>);
} // namespace o2::aod

struct MultCodes {
  static constexpr int kNoMultiplicity = 0;
  static constexpr int kMultFV0M = 1;
  static constexpr int kMultFT0M = 2;
  static constexpr int kMultFDDM = 3;
  static constexpr int kMultTracklets = 4;
  static constexpr int kMultTPC = 5;
  static constexpr int kMultNTracksPV = 6;
  static constexpr int kMultNTracksPVeta1 = 7;
  static constexpr int kCentralityFT0C = 8;
  static constexpr int kCentralityFT0M = 9;
  static constexpr int kCentralityFV0A = 10;
  static constexpr int kNMults = 10;
};

#endif // PWGLF_DATAMODEL_SPECTRATOF_H_
