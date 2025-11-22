// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ClusterHistograms.h
/// \brief Header file for histograms used in EMC cluster QA
/// \author N. Strangmann, nicolas.strangmann@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_

#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TH2.h>

#include <cstddef>
#include <string_view>

namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram
{
inline void addClusterHistograms(o2::framework::HistogramRegistry* fRegistry, bool do2DQA)
{
  fRegistry->add("Cluster/before/hE", "E_{cluster};#it{E}_{cluster} (GeV);#it{N}_{cluster}", o2::framework::kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{cluster}", o2::framework::kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", o2::framework::kTH1F, {{51, -0.5f, 50.5f}}, true);
  fRegistry->add("Cluster/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", o2::framework::kTH2F, {{280, -0.7f, 0.7f}, {180, 0, o2::constants::math::TwoPI}}, true);
  fRegistry->add("Cluster/before/hNTracks", "Number of tracks considered for TM;#it{N}_{tracks};#it{N}_{cluster}", o2::framework::kTH1F, {{20, -0.5f, 19.5}}, true);
  fRegistry->add("Cluster/before/hTrackdEtadPhi", "d#eta vs. d#varphi of matched tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {{200, -0.2f, 0.2f}, {200, -0.2f, 0.2f}}, true);

  if (do2DQA) { // Check if 2D QA histograms were selected in em-qc task
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{26, -0.5, 25.5}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, 0, 2}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{300, -150, 150}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdEta", "d#eta vs. E of matched tracks;d#eta;#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, -0.2f, 0.2f}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdPhi", "d#phi vs. E of matched tracks;d#varphi (rad.);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, -0.2f, 0.2f}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c});#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, 0., 5.}, {200, 0, 20}}, true);
  } else {
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{N}_{cluster}", o2::framework::kTH1F, {{26, -0.5, 25.5}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{N}_{cluster}", o2::framework::kTH1F, {{400, 0, 2}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{N}_{cluster}", o2::framework::kTH1F, {{600, -150, 150}}, true);
    fRegistry->add("Cluster/before/hTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c})", o2::framework::kTH1F, {{200, 0., 5.}}, true);
  }

  auto hClusterQualityCuts = fRegistry->add<TH2>("Cluster/hClusterQualityCuts", "Energy at which clusters are removed by a given cut;;#it{E} (GeV)", o2::framework::kTH2F, {{static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 2, -0.5, static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 1.5}, {500, 0, 50}}, true);
  hClusterQualityCuts->GetXaxis()->SetBinLabel(1, "In");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(2, "Definition");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(3, "Energy");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(4, "NCell");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(5, "M02");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(6, "Timing");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(7, "TM");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(8, "Sec. TM");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(9, "Exotic");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(10, "Out");

  fRegistry->addClone("Cluster/before/", "Cluster/after/");
}

template <const int cls_id, typename TCluster>
inline void fillClusterHistograms(o2::framework::HistogramRegistry* fRegistry, TCluster cluster, bool do2DQA, float weight = 1.f)
{
  static constexpr std::string_view kClusterTypes[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hE"), cluster.e(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hPt"), cluster.pt(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hEtaPhi"), cluster.eta(), cluster.phi(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hNTracks"), cluster.deltaEta().size(), weight);
  for (size_t itrack = 0; itrack < cluster.deltaEta().size(); itrack++) { // Fill TrackEtaPhi histogram with delta phi and delta eta of all tracks saved in the vectors in skimmerGammaCalo.cxx
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTrackdEtadPhi"), cluster.deltaEta()[itrack], cluster.deltaPhi()[itrack], weight);
  }
  if (do2DQA) {
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hNCell"), cluster.nCells(), cluster.e(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hM02"), cluster.m02(), cluster.e(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTime"), cluster.time(), cluster.e(), weight);
    for (size_t itrack = 0; itrack < cluster.deltaEta().size(); itrack++) {
      fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTrackEOverP"), cluster.e() / cluster.trackp()[itrack], cluster.e(), weight);
      fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTrackdEta"), cluster.deltaEta()[itrack], cluster.e(), weight);
      fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTrackdPhi"), cluster.deltaPhi()[itrack], cluster.e(), weight);
    }
  } else {
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hNCell"), cluster.nCells(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hM02"), cluster.m02(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTime"), cluster.time(), weight);
    for (size_t itrack = 0; itrack < cluster.deltaEta().size(); itrack++) {
      fRegistry->fill(HIST("Cluster/") + HIST(kClusterTypes[cls_id]) + HIST("hTrackEOverP"), cluster.e() / cluster.trackp()[itrack], weight);
    }
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram

#endif // PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
