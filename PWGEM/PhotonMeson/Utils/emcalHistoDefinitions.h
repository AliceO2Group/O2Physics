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

/// \brief commonly used histogram (axis) definitions for emcal in PWGEM
/// \author marvin.hemmer@cern.ch

#include <Framework/HistogramSpec.h>

#include <vector>

#ifndef PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_
#define PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_

using namespace o2::framework;

AxisSpec const gAxis_dEta{100, -0.2f, 0.2f, "#Delta#it{#eta}"};
AxisSpec const gAxis_dPhi{100, -0.2f, 0.2f, "#Delta#it{#varphi} (rad)"};
AxisSpec const gAxis_dR{100, 0.0f, 0.2f, "#Delta#it{R}"};
AxisSpec const gAxis_MatchedTrack{10, 0.5f, 10.5f, "matched track"};
AxisSpec const gAxis_Eta{160, -0.8f, 0.8f, "#it{#eta}"};
AxisSpec const gAxis_Phi{144, 0.f, 2 * 3.14159, "#it{#varphi} (rad)"};

std::vector<double> EClusBins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8,
                                 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8,
                                 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11,
                                 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                                 48, 49, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                                 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
AxisSpec gAxis_EClus = {EClusBins, "#it{E}_{clus} (GeV)"};

HistogramConfigSpec gHistoSpec_clusterTM_dEtadPhi({HistType::kTH2F, {gAxis_dEta, gAxis_dPhi}});
HistogramConfigSpec gHistoSpec_clusterTM_dEtaMT({HistType::kTH2F, {gAxis_dEta, {10, 0.5, 10.5}}});
HistogramConfigSpec gHistoSpec_clusterTM_dPhiMT({HistType::kTH2F, {gAxis_dPhi, {10, 0.5, 10.5}}});
HistogramConfigSpec gHistoSpec_clusterTM_dRMT({HistType::kTH2F, {gAxis_dR, {10, 0.5, 10.5}}});
HistogramConfigSpec gHistoSpec_EtaPhi({HistType::kTH2F, {gAxis_Eta, gAxis_Phi}});
HistogramConfigSpec gHistoSpec_clusterE({HistType::kTH1F, {gAxis_EClus}});
HistogramConfigSpec gHistoSpec_clusterECuts({HistType::kTH2F, {gAxis_EClus, {64, -0.5, 63.5}}});

#endif // PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_
