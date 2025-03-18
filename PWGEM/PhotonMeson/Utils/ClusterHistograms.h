// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// Header file for histograms used in EMC cluster QA
/// \author nicolas.strangmann@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_

using namespace o2::framework;

namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram
{
void addClusterHistograms(HistogramRegistry* fRegistry, bool do2DQA)
{
  fRegistry->add("Cluster/before/hE", "E_{cluster};#it{E}_{cluster} (GeV);#it{N}_{cluster}", kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{cluster}", kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", kTH1F, {{51, -0.5f, 50.5f}}, true);
  fRegistry->add("Cluster/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", kTH2F, {{280, -0.7f, 0.7f}, {180, 0, 2 * M_PI}}, true);
  fRegistry->add("Cluster/before/hNTracks", "Number of tracks considered for TM;#it{N}_{tracks};#it{N}_{cluster}", kTH1F, {{20, -0.5f, 19.5}}, true);
  fRegistry->add("Cluster/before/hTrackdEtadPhi", "d#eta vs. d#varphi of matched tracks;d#eta;d#varphi (rad.)", kTH2F, {{200, -0.2f, 0.2f}, {200, -0.2f, 0.2f}}, true);

  if (do2DQA) { // Check if 2D QA histograms were selected in em-qc task
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{E}_{cluster} (GeV)", kTH2F, {{26, -0.5, 25.5}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{E}_{cluster} (GeV)", kTH2F, {{200, 0, 2}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{E}_{cluster} (GeV)", kTH2F, {{300, -150, 150}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdEta", "d#eta vs. E of matched tracks;d#eta;#it{E}_{cluster} (GeV)", kTH2F, {{200, -0.2f, 0.2f}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdPhi", "d#phi vs. E of matched tracks;d#varphi (rad.);#it{E}_{cluster} (GeV)", kTH2F, {{200, -0.2f, 0.2f}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c});#it{E}_{cluster} (GeV)", kTH2F, {{200, 0., 5.}, {200, 0, 20}}, true);
  } else {
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{N}_{cluster}", kTH1F, {{26, -0.5, 25.5}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{N}_{cluster}", kTH1F, {{400, 0, 2}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{N}_{cluster}", kTH1F, {{600, -150, 150}}, true);
    fRegistry->add("Cluster/before/hTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c})", kTH1F, {{200, 0., 5.}}, true);
  }

  auto hClusterQualityCuts = fRegistry->add<TH2>("Cluster/hClusterQualityCuts", "Energy at which clusters are removed by a given cut;;#it{E} (GeV)", kTH2F, {{8, -0.5, 7.5}, {500, 0, 50}}, true);
  hClusterQualityCuts->GetXaxis()->SetBinLabel(1, "In");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(2, "Energy");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(3, "NCell");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(4, "M02");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(5, "Timing");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(6, "Track matching");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(7, "Exotic");
  hClusterQualityCuts->GetXaxis()->SetBinLabel(8, "Out");

  fRegistry->addClone("Cluster/before/", "Cluster/after/");
}

template <const int cls_id, typename TCluster>
void fillClusterHistograms(HistogramRegistry* fRegistry, TCluster cluster, bool do2DQA, float weight = 1.f)
{
  static constexpr std::string_view cluster_types[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hE"), cluster.e(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hPt"), cluster.pt(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hEtaPhi"), cluster.eta(), cluster.phi(), weight);
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hNTracks"), cluster.tracketa().size(), weight);
  for (size_t itrack = 0; itrack < cluster.tracketa().size(); itrack++) { // Fill TrackEtaPhi histogram with delta phi and delta eta of all tracks saved in the vectors in skimmerGammaCalo.cxx
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackdEtadPhi"), cluster.tracketa()[itrack] - cluster.eta(), cluster.trackphi()[itrack] - cluster.phi(), weight);
  }
  if (do2DQA) {
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hNCell"), cluster.nCells(), cluster.e(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hM02"), cluster.m02(), cluster.e(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTime"), cluster.time(), cluster.e(), weight);
    for (size_t itrack = 0; itrack < cluster.tracketa().size(); itrack++) {
      fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackEOverP"), cluster.e() / cluster.trackp()[itrack], cluster.e(), weight);
      fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackdEta"), cluster.tracketa()[itrack] - cluster.eta(), cluster.e(), weight);
      fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackdPhi"), cluster.trackphi()[itrack] - cluster.phi(), cluster.e(), weight);
    }
  } else {
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hNCell"), cluster.nCells(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hM02"), cluster.m02(), weight);
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTime"), cluster.time(), weight);
    for (size_t itrack = 0; itrack < cluster.tracketa().size(); itrack++) {
      fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackEOverP"), cluster.e() / cluster.trackp()[itrack], weight);
    }
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram

#endif // PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
