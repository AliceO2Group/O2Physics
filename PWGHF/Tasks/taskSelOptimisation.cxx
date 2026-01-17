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

/// \file taskSelOptimisation.cxx
/// \brief task to study preselections
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int NCutsToTestCosp = 15;
constexpr int NCutsToTestDecLen = 11;
constexpr int NCutsToTestImpParProd = 11;
constexpr int NCutsToTestMinDcAxy = 9;
constexpr int NCutsToTestMinTrackPt = 7;

constexpr float CutsCosp[NCutsToTestCosp] = {0.70, 0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995};
constexpr float CutsDecLen[NCutsToTestDecLen] = {0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1};
constexpr float CutsImpParProd[NCutsToTestImpParProd] = {-0.00005, -0.00004, -0.00003, -0.00002, -0.00001, 0., 0.00001, 0.00002, 0.00003, 0.00004, 0.00005};
constexpr float CutsMinDcAxy[NCutsToTestMinDcAxy] = {0., 0.0005, 0.001, 0.0015, 0.0020, 0.0025, 0.0030, 0.0040, 0.0050};
constexpr float CutsMinTrackPt[NCutsToTestMinTrackPt] = {0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60};

auto vecCutsCosp = std::vector<float>{CutsCosp, CutsCosp + NCutsToTestCosp};
auto vecCutsDecLen = std::vector<float>{CutsDecLen, CutsDecLen + NCutsToTestDecLen};
auto vecCutsImpParProd = std::vector<float>{CutsImpParProd, CutsImpParProd + NCutsToTestImpParProd};
auto vecCutsMinDCAxy = std::vector<float>{CutsMinDcAxy, CutsMinDcAxy + NCutsToTestMinDcAxy};
auto vecCutsMinTrackPt = std::vector<float>{CutsMinTrackPt, CutsMinTrackPt + NCutsToTestMinTrackPt};

const int n2Prong = o2::aod::hf_cand_2prong::DecayType::N2ProngDecays;
const int n3Prong = o2::aod::hf_cand_3prong::DecayType::N3ProngDecays;

constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNames2Prong = {{{"hPromptVsPtD0ToPiK", "hPromptVsPtJpsiToEE", "hPromptVsPtJpsiToMuMu", "hPromptVsPt2Prong"},
                                                                                        {"hNonPromptVsPtD0ToPiK", "hNonPromptVsPtJpsiToEE", "hNonPromptVsPtJpsiToMuMu", "hNonPromptVsPt2Prong"},
                                                                                        {"hBkgVsPtD0ToPiK", "hBkgVsPtJpsiToEE", "hBkgVsPtJpsiToMuMu", "hBkgVsPt2Prong"}}};
constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNamesCosp2Prong = {{{"hPromptCospVsPtD0ToPiK", "hPromptCospVsPtJpsiToEE", "hPromptCospVsPtJpsiToMuMu", "hPromptCospVsPt2Prong"},
                                                                                            {"hNonPromptCospVsPtD0ToPiK", "hNonPromptCospVsPtJpsiToEE", "hNonPromptCospVsPtJpsiToMuMu", "hNonPromptCospVsPt2Prong"},
                                                                                            {"hBkgCospVsPtD0ToPiK", "hBkgCospVsPtJpsiToEE", "hBkgCospVsPtJpsiToMuMu", "hBkgCospVsPt2Prong"}}};
constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNamesDecLen2Prong = {{{"hPromptDecLenVsPtD0ToPiK", "hPromptDecLenVsPtJpsiToEE", "hPromptDecLenVsPtJpsiToMuMu", "hPromptDecLenVsPt2Prong"},
                                                                                              {"hNonPromptDecLenVsPtD0ToPiK", "hNonPromptDecLenVsPtJpsiToEE", "hNonPromptDecLenVsPtJpsiToMuMu", "hNonPromptDecLenVsPt2Prong"},
                                                                                              {"hBkgDecLenVsPtD0ToPiK", "hBkgDecLenVsPtJpsiToEE", "hBkgDecLenVsPtJpsiToMuMu", "hBkgDecLenVsPt2Prong"}}};
constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNamesImpParProd2Prong = {{{"hPromptImpParProdVsPtD0ToPiK", "hPromptImpParProdVsPtJpsiToEE", "hPromptImpParProdVsPtJpsiToMuMu", "hPromptImpParProdVsPt2Prong"},
                                                                                                  {"hNonPromptImpParProdVsPtD0ToPiK", "hNonPromptImpParProdVsPtJpsiToEE", "hNonPromptImpParProdVsPtJpsiToMuMu", "hNonPromptImpParProdVsPt2Prong"},
                                                                                                  {"hBkgImpParProdVsPtD0ToPiK", "hBkgImpParProdVsPtJpsiToEE", "hBkgImpParProdVsPtJpsiToMuMu", "hBkgImpParProdVsPt2Prong"}}};
constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNamesMinDcAxy2Prong = {{{"hPromptMinDCAxyVsPtD0ToPiK", "hPromptMinDCAxyVsPtJpsiToEE", "hPromptMinDCAxyVsPtJpsiToMuMu", "hPromptMinDCAxyVsPt2Prong"},
                                                                                                {"hNonPromptMinDCAxyVsPtD0ToPiK", "hNonPromptMinDCAxyVsPtJpsiToEE", "hNonPromptMinDCAxyVsPtJpsiToMuMu", "hNonPromptMinDCAxyVsPt2Prong"},
                                                                                                {"hBkgMinDCAxyVsPtD0ToPiK", "hBkgMinDCAxyVsPtJpsiToEE", "hBkgMinDCAxyVsPtJpsiToMuMu", "hBkgMinDCAxyVsPt2Prong"}}};
constexpr std::array<std::array<std::string_view, n2Prong + 1>, 3> HistoNamesMinTrackPt2Prong = {{{"hPromptMinTrackPtVsPtD0ToPiK", "hPromptMinTrackPtVsPtJpsiToEE", "hPromptMinTrackPtVsPtJpsiToMuMu", "hPromptMinTrackPtVsPt2Prong"},
                                                                                                  {"hNonPromptMinTrackPtVsPtD0ToPiK", "hNonPromptMinTrackPtVsPtJpsiToEE", "hNonPromptMinTrackPtVsPtJpsiToMuMu", "hNonPromptMinTrackPtVsPt2Prong"},
                                                                                                  {"hBkgMinTrackPtVsPtD0ToPiK", "hBkgMinTrackPtVsPtJpsiToEE", "hBkgMinTrackPtVsPtJpsiToMuMu", "hBkgMinTrackPtVsPt2Prong"}}};

constexpr std::array<std::array<std::string_view, n3Prong + 1>, 3> HistoNames3Prong = {{{"hPromptVsPtDPlusToPiKPi", "hPromptVsPtLcToPKPi", "hPromptVsPtDsToPiKK", "hPromptVsPtXicToPKPi", "hPromptVsPt3Prong"},
                                                                                        {"hNonPromptVsPtDPlusToPiKPi", "hNonPromptVsPtLcToPKPi", "hNonPromptVsPtDsToPiKK", "hNonPromptVsPtXicToPKPi", "hNonPromptVsPt3Prong"},
                                                                                        {"hBkgVsPtDPlusToPiKPi", "hBkgVsPtLcToPKPi", "hBkgVsPtDsToPiKK", "hBkgVsPtXicToPKPi", "hBkgVsPt3Prong"}}};
constexpr std::array<std::array<std::string_view, n3Prong + 1>, 3> HistoNamesCosp3Prong = {{{"hPromptCospVsPtDPlusToPiKPi", "hPromptCospVsPtLcToPKPi", "hPromptCospVsPtDsToPiKK", "hPromptCospVsPtXicToPKPi", "hPromptCospVsPt3Prong"},
                                                                                            {"hNonPromptCospVsPtDPlusToPiKPi", "hNonPromptCospVsPtLcToPKPi", "hNonPromptCospVsPtDsToPiKK", "hNonPromptCospVsPtXicToPKPi", "hNonPromptCospVsPt3Prong"},
                                                                                            {"hBkgCospVsPtDPlusToPiKPi", "hBkgCospVsPtLcToPKPi", "hBkgCospVsPtDsToPiKK", "hBkgCospVsPtXicToPKPi", "hBkgCospVsPt3Prong"}}};
constexpr std::array<std::array<std::string_view, n3Prong + 1>, 3> HistoNamesDecLen3Prong = {{{"hPromptDecLenVsPtDPlusToPiKPi", "hPromptDecLenVsPtLcToPKPi", "hPromptDecLenVsPtDsToPiKK", "hPromptDecLenVsPtXicToPKPi", "hPromptDecLenVsPt3Prong"},
                                                                                              {"hNonPromptDecLenVsPtDPlusToPiKPi", "hNonPromptDecLenVsPtLcToPKPi", "hNonPromptDecLenVsPtDsToPiKK", "hNonPromptDecLenVsPtXicToPKPi", "hNonPromptDecLenVsPt3Prong"},
                                                                                              {"hBkgDecLenVsPtDPlusToPiKPi", "hBkgDecLenVsPtLcToPKPi", "hBkgDecLenVsPtDsToPiKK", "hBkgDecLenVsPtXicToPKPi", "hBkgDecLenVsPt3Prong"}}};
constexpr std::array<std::array<std::string_view, n3Prong + 1>, 3> HistoNamesMinDcAxy3Prong = {{{"hPromptMinDCAxyVsPtDPlusToPiKPi", "hPromptMinDCAxyVsPtLcToPKPi", "hPromptMinDCAxyVsPtDsToPiKK", "hPromptMinDCAxyVsPtXicToPKPi", "hPromptMinDCAxyVsPt3Prong"},
                                                                                                {"hNonPromptMinDCAxyVsPtDPlusToPiKPi", "hNonPromptMinDCAxyVsPtLcToPKPi", "hNonPromptMinDCAxyVsPtDsToPiKK", "hNonPromptMinDCAxyVsPtXicToPKPi", "hNonPromptMinDCAxyVsPt3Prong"},
                                                                                                {"hBkgMinDCAxyVsPtDPlusToPiKPi", "hBkgMinDCAxyVsPtLcToPKPi", "hBkgMinDCAxyVsPtDsToPiKK", "hBkgMinDCAxyVsPtXicToPKPi", "hBkgMinDCAxyVsPt3Prong"}}};
constexpr std::array<std::array<std::string_view, n3Prong + 1>, 3> HistoNamesMinTrackPt3Prong = {{{"hPromptMinTrackPtVsPtDPlusToPiKPi", "hPromptMinTrackPtVsPtLcToPKPi", "hPromptMinTrackPtVsPtDsToPiKK", "hPromptMinTrackPtVsPtXicToPKPi", "hPromptMinTrackPtVsPt3Prong"},
                                                                                                  {"hNonPromptMinTrackPtVsPtDPlusToPiKPi", "hNonPromptMinTrackPtVsPtLcToPKPi", "hNonPromptMinTrackPtVsPtDsToPiKK", "hNonPromptMinTrackPtVsPtXicToPKPi", "hNonPromptMinTrackPtVsPt3Prong"},
                                                                                                  {"hBkgMinTrackPtVsPtDPlusToPiKPi", "hBkgMinTrackPtVsPtLcToPKPi", "hBkgMinTrackPtVsPtDsToPiKK", "hBkgMinTrackPtVsPtXicToPKPi", "hBkgMinTrackPtVsPt3Prong"}}};

std::array<std::array<std::shared_ptr<TH1>, n2Prong + 1>, 3> histPt2Prong{};
std::array<std::array<std::shared_ptr<TH2>, n2Prong + 1>, 3> histCospVsPt2Prong{};
std::array<std::array<std::shared_ptr<TH2>, n2Prong + 1>, 3> histDecLenVsPt2Prong{};
std::array<std::array<std::shared_ptr<TH2>, n2Prong + 1>, 3> histImpParProdVsPt2Prong{};
std::array<std::array<std::shared_ptr<TH2>, n2Prong + 1>, 3> histMinDCAxyVsPt2Prong{};
std::array<std::array<std::shared_ptr<TH2>, n2Prong + 1>, 3> histMinTrackPtVsPt2Prong{};

std::array<std::array<std::shared_ptr<TH1>, n3Prong + 1>, 3> histPt3Prong{};
std::array<std::array<std::shared_ptr<TH2>, n3Prong + 1>, 3> histCospVsPt3Prong{};
std::array<std::array<std::shared_ptr<TH2>, n3Prong + 1>, 3> histDecLenVsPt3Prong{};
std::array<std::array<std::shared_ptr<TH2>, n3Prong + 1>, 3> histMinDCAxyVsPt3Prong{};
std::array<std::array<std::shared_ptr<TH2>, n3Prong + 1>, 3> histMinTrackPtVsPt3Prong{};

} // namespace

struct HfSelOptimisation {
  Configurable<std::vector<float>> cutsToTestCpa{"cutsToTestCpa", std::vector<float>{vecCutsCosp}, "cos(theta_P) cut values to test"};
  Configurable<std::vector<float>> cutsToTestDecLen{"cutsToTestDecLen", std::vector<float>{vecCutsDecLen}, "decay length cut values to test"};
  Configurable<std::vector<float>> cutsToTestImpParProd{"cutsToTestImpParProd", std::vector<float>{vecCutsImpParProd}, "impact parameter product cut values to test (2-prongs only)"};
  Configurable<std::vector<float>> cutsToTestMinDcaXY{"cutsToTestMinDcaXY", std::vector<float>{vecCutsMinDCAxy}, "min DCA xy cut values to test"};
  Configurable<std::vector<float>> cutsToTestMinTrackPt{"cutsToTestMinTrackPt", std::vector<float>{vecCutsMinTrackPt}, "min track pT cut values to test"};

  ConfigurableAxis ptBinning{"ptBinning", {0, 0., 2., 5., 20.}, "pT bin limits"};

  AxisSpec axisPt = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
  // quantized axes
  AxisSpec axisCosp = {static_cast<int>(cutsToTestCpa->size()), 0.5, cutsToTestCpa->size() + 0.5, "cos(#theta_{P}) >"};
  AxisSpec axisDecLen = {static_cast<int>(cutsToTestDecLen->size()), 0.5, cutsToTestDecLen->size() + 0.5, "decay length (cm) >"};
  AxisSpec axisImpParProd = {static_cast<int>(cutsToTestImpParProd->size()), 0.5, cutsToTestImpParProd->size() + 0.5, "#it{d}_{0}#times#it{d}_{0} (cm^{2}) <"};
  AxisSpec axisMinDCAxy = {static_cast<int>(cutsToTestMinDcaXY->size()), 0.5, cutsToTestMinDcaXY->size() + 0.5, "min track #it{d}_{0} (cm) >"};
  AxisSpec axisMinTrackPt = {static_cast<int>(cutsToTestMinTrackPt->size()), 0.5, cutsToTestMinTrackPt->size() + 0.5, "min track #it{p}_{T} (cm) >"};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext const&)
  {
    for (int iOrig{0}; iOrig < 3; iOrig++) {
      for (int i2Prong = 0; i2Prong < n2Prong + 1; ++i2Prong) {
        histPt2Prong[iOrig][i2Prong] = registry.add<TH1>(HistoNames2Prong[iOrig][i2Prong].data(), "", HistType::kTH1F, {axisPt});
        histCospVsPt2Prong[iOrig][i2Prong] = registry.add<TH2>(HistoNamesCosp2Prong[iOrig][i2Prong].data(), "", HistType::kTH2F, {axisPt, axisCosp});
        for (int iBin{0}; iBin < histCospVsPt2Prong[iOrig][i2Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histCospVsPt2Prong[iOrig][i2Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestCpa->at(iBin)));
        }
        histDecLenVsPt2Prong[iOrig][i2Prong] = registry.add<TH2>(HistoNamesDecLen2Prong[iOrig][i2Prong].data(), "", HistType::kTH2F, {axisPt, axisDecLen});
        for (int iBin{0}; iBin < histDecLenVsPt2Prong[iOrig][i2Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histDecLenVsPt2Prong[iOrig][i2Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.3f", cutsToTestDecLen->at(iBin)));
        }
        histImpParProdVsPt2Prong[iOrig][i2Prong] = registry.add<TH2>(HistoNamesImpParProd2Prong[iOrig][i2Prong].data(), "", HistType::kTH2F, {axisPt, axisImpParProd});
        for (int iBin{0}; iBin < histImpParProdVsPt2Prong[iOrig][i2Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histImpParProdVsPt2Prong[iOrig][i2Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestImpParProd->at(iBin)));
        }
        histMinDCAxyVsPt2Prong[iOrig][i2Prong] = registry.add<TH2>(HistoNamesMinDcAxy2Prong[iOrig][i2Prong].data(), "", HistType::kTH2F, {axisPt, axisMinDCAxy});
        for (int iBin{0}; iBin < histMinDCAxyVsPt2Prong[iOrig][i2Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histMinDCAxyVsPt2Prong[iOrig][i2Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestMinDcaXY->at(iBin)));
        }
        histMinTrackPtVsPt2Prong[iOrig][i2Prong] = registry.add<TH2>(HistoNamesMinTrackPt2Prong[iOrig][i2Prong].data(), "", HistType::kTH2F, {axisPt, axisMinTrackPt});
        for (int iBin{0}; iBin < histMinTrackPtVsPt2Prong[iOrig][i2Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histMinTrackPtVsPt2Prong[iOrig][i2Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.2f", cutsToTestMinTrackPt->at(iBin)));
        }
      }
      for (int i3Prong{0}; i3Prong < n3Prong + 1; ++i3Prong) {
        histPt3Prong[iOrig][i3Prong] = registry.add<TH1>(HistoNames3Prong[iOrig][i3Prong].data(), "", HistType::kTH1F, {axisPt});
        histCospVsPt3Prong[iOrig][i3Prong] = registry.add<TH2>(HistoNamesCosp3Prong[iOrig][i3Prong].data(), "", HistType::kTH2F, {axisPt, axisCosp});
        for (int iBin{0}; iBin < histCospVsPt3Prong[iOrig][i3Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histCospVsPt3Prong[iOrig][i3Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestCpa->at(iBin)));
        }
        histDecLenVsPt3Prong[iOrig][i3Prong] = registry.add<TH2>(HistoNamesDecLen3Prong[iOrig][i3Prong].data(), "", HistType::kTH2F, {axisPt, axisDecLen});
        for (int iBin{0}; iBin < histDecLenVsPt3Prong[iOrig][i3Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histDecLenVsPt3Prong[iOrig][i3Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestDecLen->at(iBin)));
        }
        histMinDCAxyVsPt3Prong[iOrig][i3Prong] = registry.add<TH2>(HistoNamesMinDcAxy3Prong[iOrig][i3Prong].data(), "", HistType::kTH2F, {axisPt, axisMinDCAxy});
        for (int iBin{0}; iBin < histMinDCAxyVsPt3Prong[iOrig][i3Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histMinDCAxyVsPt3Prong[iOrig][i3Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestMinDcaXY->at(iBin)));
        }
        histMinTrackPtVsPt3Prong[iOrig][i3Prong] = registry.add<TH2>(HistoNamesMinTrackPt3Prong[iOrig][i3Prong].data(), "", HistType::kTH2F, {axisPt, axisMinTrackPt});
        for (int iBin{0}; iBin < histMinTrackPtVsPt3Prong[iOrig][i3Prong]->GetYaxis()->GetNbins(); ++iBin) {
          histMinTrackPtVsPt3Prong[iOrig][i3Prong]->GetYaxis()->SetBinLabel(iBin + 1, Form("%0.4f", cutsToTestMinTrackPt->at(iBin)));
        }
      }
    }
  }

  /// Conjugate-dependent topological cuts
  /// \param candType is the candidate channel
  /// \param candOrig is candidate type (Prompt, NonPrompt, Bkg)
  /// \param candidate is a candidate
  /// \param tracks is the array of daughter tracks
  template <std::size_t CandType, std::size_t CandOrig, typename T1, typename T2>
  void testSelections2Prong(const T1& candidate, const T2& tracks)
  {
    auto pT = candidate.pt();
    std::array<double, 2> absDCA{std::abs(tracks[0].dcaXY()), std::abs(tracks[1].dcaXY())};
    std::sort(absDCA.begin(), absDCA.end());

    std::array<double, 2> ptTrack{tracks[0].pt(), tracks[1].pt()};
    std::sort(ptTrack.begin(), ptTrack.end());

    histPt2Prong[CandOrig][CandType]->Fill(pT);

    for (std::size_t iCospCut{0}; iCospCut < cutsToTestCpa->size(); ++iCospCut) {
      if (candidate.cpa() > cutsToTestCpa->at(iCospCut)) {
        histCospVsPt2Prong[CandOrig][CandType]->Fill(pT, iCospCut + 1);
      }
    }

    for (std::size_t iDecLenCut{0}; iDecLenCut < cutsToTestDecLen->size(); ++iDecLenCut) {
      if (candidate.decayLength() > cutsToTestDecLen->at(iDecLenCut)) {
        histDecLenVsPt2Prong[CandOrig][CandType]->Fill(pT, iDecLenCut + 1);
      }
    }

    for (std::size_t iImpParProd{0}; iImpParProd < cutsToTestImpParProd->size(); ++iImpParProd) {
      if (candidate.impactParameterProduct() < cutsToTestImpParProd->at(iImpParProd)) {
        histImpParProdVsPt2Prong[CandOrig][CandType]->Fill(pT, iImpParProd + 1);
      }
    }

    for (std::size_t iMinDCAxy{0}; iMinDCAxy < cutsToTestMinDcaXY->size(); ++iMinDCAxy) {
      if (absDCA[0] > cutsToTestMinDcaXY->at(iMinDCAxy)) {
        histMinDCAxyVsPt2Prong[CandOrig][CandType]->Fill(pT, iMinDCAxy + 1);
      }
    }

    for (std::size_t iMinTrackPt{0}; iMinTrackPt < cutsToTestMinTrackPt->size(); ++iMinTrackPt) {
      if (ptTrack[0] > cutsToTestMinTrackPt->at(iMinTrackPt)) {
        histMinTrackPtVsPt2Prong[CandOrig][CandType]->Fill(pT, iMinTrackPt + 1);
      }
    }
  }

  /// Conjugate-dependent topological cuts
  /// \param candType is the candidate channel
  /// \param candOrig is candidate type (Prompt, NonPrompt, Bkg)
  /// \param candidate is a candidate
  /// \param tracks is the array of doughter tracks
  template <std::size_t CandType, std::size_t CandOrig, typename T1, typename T2>
  void testSelections3Prong(const T1& candidate, const T2& tracks)
  {
    auto pT = candidate.pt();
    std::array<double, 3> absDCA{std::abs(tracks[0].dcaXY()), std::abs(tracks[1].dcaXY()), std::abs(tracks[2].dcaXY())};
    std::sort(absDCA.begin(), absDCA.end());

    std::array<double, 3> ptTrack{tracks[0].pt(), tracks[1].pt(), tracks[2].pt()};
    std::sort(ptTrack.begin(), ptTrack.end());

    histPt3Prong[CandOrig][CandType]->Fill(pT);

    for (std::size_t iCospCut{0}; iCospCut < cutsToTestCpa->size(); ++iCospCut) {
      if (candidate.cpa() > cutsToTestCpa->at(iCospCut)) {
        histCospVsPt3Prong[CandOrig][CandType]->Fill(pT, iCospCut + 1);
      }
    }

    for (std::size_t iDecLenCut{0}; iDecLenCut < cutsToTestDecLen->size(); ++iDecLenCut) {
      if (candidate.decayLength() > cutsToTestDecLen->at(iDecLenCut)) {
        histDecLenVsPt3Prong[CandOrig][CandType]->Fill(pT, iDecLenCut + 1);
      }
    }

    for (std::size_t iMinDCAxy{0}; iMinDCAxy < cutsToTestMinDcaXY->size(); ++iMinDCAxy) {
      if (absDCA[0] > cutsToTestMinDcaXY->at(iMinDCAxy)) {
        histMinDCAxyVsPt3Prong[CandOrig][CandType]->Fill(pT, iMinDCAxy + 1);
      }
    }

    for (std::size_t iMinTrackPt{0}; iMinTrackPt < cutsToTestMinTrackPt->size(); ++iMinTrackPt) {
      if (ptTrack[0] > cutsToTestMinTrackPt->at(iMinTrackPt)) {
        histMinTrackPtVsPt3Prong[CandOrig][CandType]->Fill(pT, iMinTrackPt + 1);
      }
    }
  }

  void process(soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec> const& cand2Prongs,
               soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec> const& cand3Prongs,
               aod::TracksWDca const&)
  {
    // looping over 2-prong candidates
    for (const auto& cand2Prong : cand2Prongs) {

      auto trackPos = cand2Prong.prong0_as<aod::TracksWDca>(); // positive daughter
      auto trackNeg = cand2Prong.prong1_as<aod::TracksWDca>(); // negative daughter
      std::array const tracks = {trackPos, trackNeg};

      bool isPrompt = false, isNonPrompt = false, isBkg = false;
      for (int iDecay{0}; iDecay < n2Prong; ++iDecay) {
        if (TESTBIT(cand2Prong.hfflag(), iDecay)) {
          if (std::abs(cand2Prong.flagMcMatchRec()) == BIT(iDecay)) { // FIXME: Migrate to DecayChannelMain
            if (cand2Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
              isPrompt = true;
              switch (iDecay) {
                case o2::aod::hf_cand_2prong::DecayType::D0ToPiK:
                  testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::D0ToPiK, 0>(cand2Prong, tracks);
                  break;
                case o2::aod::hf_cand_2prong::DecayType::JpsiToEE:
                  testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::JpsiToEE, 0>(cand2Prong, tracks);
                  break;
              }
            } else if (cand2Prong.originMcRec() == RecoDecay::OriginType::NonPrompt) {
              isNonPrompt = true;
              switch (iDecay) {
                case o2::aod::hf_cand_2prong::DecayType::D0ToPiK:
                  testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::D0ToPiK, 1>(cand2Prong, tracks);
                  break;
                case o2::aod::hf_cand_2prong::DecayType::JpsiToEE:
                  testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::JpsiToEE, 1>(cand2Prong, tracks);
                  break;
              }
            }
          } else {
            isBkg = true;
            switch (iDecay) {
              case o2::aod::hf_cand_2prong::DecayType::D0ToPiK:
                testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::D0ToPiK, 2>(cand2Prong, tracks);
                break;
              case o2::aod::hf_cand_2prong::DecayType::JpsiToEE:
                testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::JpsiToEE, 2>(cand2Prong, tracks);
                break;
            }
          }
        }
      }

      if (isPrompt) {
        testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::N2ProngDecays, 0>(cand2Prong, tracks);
      } else if (isNonPrompt) {
        testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::N2ProngDecays, 1>(cand2Prong, tracks);
      } else if (isBkg) {
        testSelections2Prong<o2::aod::hf_cand_2prong::DecayType::N2ProngDecays, 2>(cand2Prong, tracks);
      }
    } // loop over 2-prong candidates

    // looping over 3-prong candidates
    for (const auto& cand3Prong : cand3Prongs) {

      auto trackFirst = cand3Prong.prong0_as<aod::TracksWDca>();  // first daughter
      auto trackSecond = cand3Prong.prong1_as<aod::TracksWDca>(); // second daughter
      auto trackThird = cand3Prong.prong2_as<aod::TracksWDca>();  // third daughter
      std::array const tracks = {trackFirst, trackSecond, trackThird};

      bool isPrompt = false, isNonPrompt = false, isBkg = false;
      for (int iDecay{0}; iDecay < n3Prong; ++iDecay) {
        if (TESTBIT(cand3Prong.hfflag(), iDecay)) {
          if (std::abs(cand3Prong.flagMcMatchRec()) == BIT(iDecay)) { // FIXME: Migrate to DecayChannelMain
            if (cand3Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
              isPrompt = true;
              switch (iDecay) {
                case o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi, 0>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::LcToPKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::LcToPKPi, 0>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::DsToKKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DsToKKPi, 0>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::XicToPKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::XicToPKPi, 0>(cand3Prong, tracks);
                  break;
              }
            } else if (cand3Prong.originMcRec() == RecoDecay::OriginType::NonPrompt) {
              isNonPrompt = true;
              switch (iDecay) {
                case o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi, 1>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::LcToPKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::LcToPKPi, 1>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::DsToKKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DsToKKPi, 1>(cand3Prong, tracks);
                  break;
                case o2::aod::hf_cand_3prong::DecayType::XicToPKPi:
                  testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::XicToPKPi, 1>(cand3Prong, tracks);
                  break;
              }
            }
          } else {
            isBkg = true;
            switch (iDecay) {
              case o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi:
                testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi, 2>(cand3Prong, tracks);
                break;
              case o2::aod::hf_cand_3prong::DecayType::LcToPKPi:
                testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::LcToPKPi, 2>(cand3Prong, tracks);
                break;
              case o2::aod::hf_cand_3prong::DecayType::DsToKKPi:
                testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::DsToKKPi, 2>(cand3Prong, tracks);
                break;
              case o2::aod::hf_cand_3prong::DecayType::XicToPKPi:
                testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::XicToPKPi, 2>(cand3Prong, tracks);
                break;
            }
          }
        }
      }
      if (isPrompt) {
        testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::N3ProngDecays, 0>(cand3Prong, tracks);
      } else if (isNonPrompt) {
        testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::N3ProngDecays, 1>(cand3Prong, tracks);
      } else if (isBkg) {
        testSelections3Prong<o2::aod::hf_cand_3prong::DecayType::N3ProngDecays, 2>(cand3Prong, tracks);
      }
    } // loop over 3-prong candidates
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfSelOptimisation>(cfgc)};
}
