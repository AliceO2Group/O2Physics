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
#include <vector>
#include <memory>
#include "CommonConstants/MathConstants.h"

using namespace o2;

JFFlucAnalysisO2Hist::JFFlucAnalysisO2Hist(HistogramRegistry& registry, AxisSpec& axisMultiplicity, AxisSpec& phiAxis, AxisSpec& etaAxis, AxisSpec& zvtAxis, AxisSpec& ptAxis, AxisSpec& massAxis, const TString& folder) : JFFlucAnalysis()
{

  ph1[HIST_TH1_CENTRALITY] = std::get<std::shared_ptr<TH1>>(registry.add(Form("%s/h_cent", folder.Data()), "multiplicity/centrality", {HistType::kTH1F, {axisMultiplicity}})).get();
  ph1[HIST_TH1_IMPACTPARAM] = std::get<std::shared_ptr<TH1>>(registry.add(Form("%s/h_IP", folder.Data()), "impact parameter", {HistType::kTH1F, {{400, -2.0, 20.0}}})).get();
  ph1[HIST_TH1_ZVERTEX] = std::get<std::shared_ptr<TH1>>(registry.add(Form("%s/h_vertex", folder.Data()), "z vertex", {HistType::kTH1F, {{100, -20.0, 20.0}}})).get();

  AxisSpec chgAxis = {3, -1.5, 1.5, "charge"};
  AxisSpec typeAxis = {2, -0.5, 1.5, "type"};
  pht[HIST_THN_PHIETAZ] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_phietaz", folder.Data()), "multiplicity/centrality, type, phi, eta, z", {HistType::kTHnSparseF, {axisMultiplicity, typeAxis, phiAxis, etaAxis, zvtAxis}})).get();
  pht[HIST_THN_PTETA] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_pteta", folder.Data()), "(corrected) multiplicity/centrality, pT, eta, charge", {HistType::kTHnSparseF, {axisMultiplicity, ptAxis, etaAxis, chgAxis}})).get();
  pht[HIST_THN_PHIETA] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_phieta", folder.Data()), "(corrected) multiplicity/centrality, phi, eta", {HistType::kTHnSparseF, {axisMultiplicity, phiAxis, etaAxis}})).get();
  AxisSpec hAxis = {kNH, -0.5, static_cast<double>(kNH - 1) + 0.5, "#it{n}"};
  AxisSpec kAxis = {nKL, -0.5, static_cast<double>(nKL - 1) + 0.5, "#it{k}"};
  AxisSpec vnAxis = {2048, -0.1, 0.1, "#it{V}_#it{n}"};
  pht[HIST_THN_SC_with_QC_4corr] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_SC_with_QC_4corr", folder.Data()), "SC_with_QC_4corr", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, hAxis, hAxis, {2048, -0.001, 0.001, "correlation"}}})).get();
  pht[HIST_THN_SC_with_QC_2corr] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_SC_with_QC_2corr", folder.Data()), "SC_with_QC_2corr", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, hAxis, {2048, -0.1, 0.1, "correlation"}}})).get();
  pht[HIST_THN_SC_with_QC_2corr_gap] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_SC_with_QC_2corr_gap", folder.Data()), "SC_with_QC_2corr_gap", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, hAxis, {2048, -0.1, 0.1, "correlation"}}})).get();
  for (UInt_t i = HIST_THN_V4V2star_2; i < HIST_THN_COUNT; ++i)
    pht[i] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/h_corrC%02u", folder.Data(), i - HIST_THN_V4V2star_2), "correlator", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, {2048, -3.0, 3.0, "correlation"}}})).get();
  for (UInt_t i = 0; i < HIST_THN_COUNT; ++i)
    pht[i]->Sumw2();

  phs[HIST_THN_SPARSE_VN] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/hvna", folder.Data()), "#it{V}_#it{n}^#it{k}", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, hAxis, kAxis, vnAxis}})).get();
  phs[HIST_THN_SPARSE_VN_VN] = std::get<std::shared_ptr<THnSparse>>(registry.add(Form("%s/hvn_vn", folder.Data()), "#it{V}_#it{n_1}^#it{k_1}#it{V}_#it{n_2}^#it{k_2}", {HistType::kTHnSparseF, {axisMultiplicity, massAxis, hAxis, kAxis, hAxis, kAxis, vnAxis}})).get();
  for (UInt_t i = 0; i < HIST_THN_SPARSE_COUNT; ++i)
    phs[i]->Sumw2();
}

JFFlucAnalysisO2Hist::~JFFlucAnalysisO2Hist()
{
  //
}
