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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \since Sep 2024

#include "JFFlucAnalysisO2Hist.h"

using namespace o2;
using namespace o2::framework;
// using namespace o2::framework::expressions;

JFFlucAnalysisO2Hist::JFFlucAnalysisO2Hist(HistogramRegistry& registry, AxisSpec& axisMultiplicity) : JFFlucAnalysis()
{
  //
  std::vector<double> ptBinning = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0};
  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
  AxisSpec phiAxis = {50, -TMath::Pi(), TMath::Pi(), "#phi"};
  AxisSpec etaAxis = {40, -2.0, 2.0, "#eta"};
  AxisSpec zvtAxis = {20, -10.0, 10.0, "zvtx"};
  pht[HIST_THN_PHIETAZ] = std::get<std::shared_ptr<THn>>(registry.add("h_phietaz", "multiplicity/centrality, phi, eta, z", {HistType::kTHnF, {axisMultiplicity, phiAxis, etaAxis, zvtAxis}})).get();
  pht[HIST_THN_PTETA] = std::get<std::shared_ptr<THn>>(registry.add("h_pteta", "(corrected) multiplicity/centrality, pT, eta", {HistType::kTHnF, {axisMultiplicity, ptAxis, etaAxis}})).get();
  pht[HIST_THN_PHIETA] = std::get<std::shared_ptr<THn>>(registry.add("h_phieta", "(corrected) multiplicity/centrality, phi, eta", {HistType::kTHnF, {axisMultiplicity, phiAxis, etaAxis}})).get();
  AxisSpec hAxis = {kNH, -0.5, (double)kNH + 0.5, "#it{n}"};
  AxisSpec kAxis = {nKL, -0.5, (double)nKL + 0.5, "#it{k}"};
  AxisSpec vnAxis = {1024, -1.0, 1.0, "#it{V}_#it{n}"};
  pht[HIST_THN_VN] = std::get<std::shared_ptr<THn>>(registry.add("hvna", "#it{V}_#it{n}^#it{k}", {HistType::kTHnF, {axisMultiplicity, hAxis, kAxis, vnAxis}})).get();
  pht[HIST_THN_VN_VN] = std::get<std::shared_ptr<THn>>(registry.add("hvn_vn", "#it{V}_#it{n_1}^#it{k_1}#it{V}_#it{n_2}^#it{k_2}", {HistType::kTHnF, {axisMultiplicity, hAxis, kAxis, hAxis, kAxis, vnAxis}})).get();
  for (UInt_t i = 0; i < HIST_THN_COUNT; ++i)
    pht[i]->Sumw2();

  ph1[HIST_TH1_CENTRALITY] = std::get<std::shared_ptr<TH1>>(registry.add("h_cent", "multiplicity/centrality", {HistType::kTH1F, {axisMultiplicity}})).get();
  ph1[HIST_TH1_IMPACTPARAM] = std::get<std::shared_ptr<TH1>>(registry.add("h_IP", "impact parameter", {HistType::kTH1F, {{400, -2.0, 20.0}}})).get();
  ph1[HIST_TH1_ZVERTEX] = std::get<std::shared_ptr<TH1>>(registry.add("h_vertex", "z vertex", {HistType::kTH1F, {{100, -20.0, 20.0}}})).get();
}

JFFlucAnalysisO2Hist::~JFFlucAnalysisO2Hist()
{
  //
}
