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

/// \file emcalHistoDefinitions.h
/// \brief commonly used histograms and axes definitions for emcal in PWGEM
/// \author marvin.hemmer@cern.ch

#include <Framework/HistogramSpec.h>

#include <vector>

#ifndef PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_
#define PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_

o2::framework::AxisSpec const gAxisdEta{200, -0.1f, 0.1f, "#Delta#it{#eta}"};
o2::framework::AxisSpec const gAxisdPhi{200, -0.1f, 0.1f, "#Delta#it{#varphi} (rad)"};
o2::framework::AxisSpec const gAxisdR{100, 0.0f, 0.2f, "#Delta#it{R}"};
o2::framework::AxisSpec const gAxisMatchedTrack{10, 0.5f, 10.5f, "matched track"};
o2::framework::AxisSpec const gAxisEta{160, -0.8f, 0.8f, "#it{#eta}"};
o2::framework::AxisSpec const gAxisPhi{144, 0.f, 2 * 3.14159, "#it{#varphi} (rad)"};
o2::framework::AxisSpec const gAxisEoverP{250, 0.f, 2.5f, "#it{E}_{clus.}/#it{p}_{track} (#it{c})"};

const std::vector<double> eClusBins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                       1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8,
                                       3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8,
                                       5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11,
                                       12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                       24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                       36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                                       48, 49, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                                       100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
o2::framework::AxisSpec const gAxisEClus = {eClusBins, "#it{E}_{clus} (GeV)"};

o2::framework::HistogramConfigSpec const gHistoSpecClusterTMdEtadPhi({o2::framework::HistType::kTH2F, {gAxisdEta, gAxisdPhi}});
o2::framework::HistogramConfigSpec const gHistoSpecClusterTMdEtaMT({o2::framework::HistType::kTH2F, {gAxisdEta, {10, 0.5, 10.5}}});
o2::framework::HistogramConfigSpec const gHistoSpecClusterTMdPhiMT({o2::framework::HistType::kTH2F, {gAxisdPhi, {10, 0.5, 10.5}}});
o2::framework::HistogramConfigSpec const gHistoSpecClusterTMdRMT({o2::framework::HistType::kTH2F, {gAxisdR, {10, 0.5, 10.5}}});
o2::framework::HistogramConfigSpec const gHistoSpecEtaPhi({o2::framework::HistType::kTH2F, {gAxisEta, gAxisPhi}});
o2::framework::HistogramConfigSpec const gHistoSpecClusterE({o2::framework::HistType::kTH1F, {gAxisEClus}});
o2::framework::HistogramConfigSpec const gHistoSpecClusterECuts({o2::framework::HistType::kTH2F, {gAxisEClus, {64, -0.5, 63.5}}});
o2::framework::HistogramConfigSpec const gHistoSpecTMEoverP({o2::framework::HistType::kTH2D, {gAxisEClus, gAxisEoverP}});

#endif // PWGEM_PHOTONMESON_UTILS_EMCALHISTODEFINITIONS_H_
