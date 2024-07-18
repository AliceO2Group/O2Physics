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
  fRegistry->add("Cluster/before/hE", "E_{cluster};#it{E}_{cluster} (GeV);#it{N}_{cluster}", kTH1F, {{500, 0, 50}}, true);
  fRegistry->add("Cluster/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{cluster}", kTH1F, {{1000, 0.0f, 20}}, true);
  fRegistry->add("Cluster/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", kTH1F, {{101, -0.5f, 100.5f}}, true);
  fRegistry->add("Cluster/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", kTH2F, {{200, -1.0f, 1.0f}, {180, 0, 2 * M_PI}}, true);
  fRegistry->add("Cluster/before/hTrackEtaPhi", "d#eta vs. d#varphi of matched tracks;d#eta;d#varphi (rad.)", kTH2F, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}, true);

  if (do2DQA) { // Check if 2D QA histograms were selected in em-qc task
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{E}_{cluster} (GeV)", kTH2F, {{51, -0.5, 50.5}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{E}_{cluster} (GeV)", kTH2F, {{500, 0, 5}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{E}_{cluster} (GeV)", kTH2F, {{100, -250, 250}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hCellTime", "Cell time;#it{t}_{cell} (ns);#it{E}_{cluster} (GeV)", kTH2F, {{100, -250, 250}, {200, 0, 20}}, true);
  } else {
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{N}_{cluster}", kTH1F, {{51, -0.5, 50.5}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{N}_{cluster}", kTH1F, {{500, 0, 5}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{N}_{cluster}", kTH1F, {{500, -250, 250}}, true);
    fRegistry->add("Cluster/before/hCellTime", "Cluster time;#it{t}_{cell} (ns);#it{N}_{cluster}", kTH1F, {{500, -250, 250}}, true);
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

template <const int cls_id>
void fillClusterHistograms(HistogramRegistry* fRegistry, SkimEMCCluster cluster, bool do2DQA)
{
  static constexpr std::string_view cluster_types[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hE"), cluster.pt());
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hPt"), cluster.e());
  fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hEtaPhi"), cluster.eta(), cluster.phi());
  for (size_t itrack = 0; itrack < cluster.tracketa().size(); itrack++) { // Fill TrackEtaPhi histogram with delta phi and delta eta of all tracks saved in the vectors in skimmerGammaCalo.cxx
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTrackEtaPhi"), cluster.tracketa()[itrack] - cluster.eta(), cluster.trackphi()[itrack] - cluster.phi());
  }
  if (do2DQA) {
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hNCell"), cluster.nCells(), cluster.e());
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hM02"), cluster.m02(), cluster.e());
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTime"), cluster.time(), cluster.e());
  } else {
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hNCell"), cluster.nCells());
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hM02"), cluster.m02());
    fRegistry->fill(HIST("Cluster/") + HIST(cluster_types[cls_id]) + HIST("hTime"), cluster.time());
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram

#endif // PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
