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
#include <Framework/ASoA.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TH2.h>

#include <cstddef>

namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram
{
inline void addClusterHistograms(o2::framework::HistogramRegistry* fRegistry, bool do2DQA)
{
  const o2::framework::AxisSpec thAxisMomentum{250, 0., 25., "#it{p}_{T} (GeV/#it{c})"};
  const o2::framework::AxisSpec thAxisDEta{200, -0.1, 0.1, "#Delta#eta"};
  const o2::framework::AxisSpec thAxisDPhi{200, -0.1, 0.1, "#Delta#varphi (rad)"};
  const o2::framework::AxisSpec thAxisEnergy{500, 0., 50., "#it{E} (GeV)"};
  const o2::framework::AxisSpec thAxisEta{320, -0.8, 0.8, "#eta"};
  const o2::framework::AxisSpec thAxisPhi{500, 0, o2::constants::math::TwoPI, "#varphi (rad)"};

  fRegistry->add("Cluster/before/hE", "E_{cluster};#it{E}_{cluster} (GeV);#it{N}_{cluster}", o2::framework::kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hPt", "Transverse momenta of clusters;#it{p}_{T} (GeV/c);#it{N}_{cluster}", o2::framework::kTH1F, {{500, 0.0f, 50}}, true);
  fRegistry->add("Cluster/before/hNgamma", "Number of #gamma candidates per collision;#it{N}_{#gamma} per collision;#it{N}_{collisions}", o2::framework::kTH1F, {{1001, -0.5f, 1000.5f}}, true);
  fRegistry->add("Cluster/before/hEtaPhi", "#eta vs #varphi;#eta;#varphi (rad.)", o2::framework::kTH2F, {thAxisEta, thAxisPhi}, true);
  fRegistry->add("Cluster/before/hNTracks", "Number of tracks considered for TM;#it{N}_{tracks};#it{N}_{cluster}", o2::framework::kTH1F, {{20, -0.5f, 19.5}}, true);
  fRegistry->add("Cluster/before/hNSecTracks", "Number of secondary tracks considered for TM;#it{N}_{tracks};#it{N}_{cluster}", o2::framework::kTH1F, {{20, -0.5f, 19.5}}, true);
  fRegistry->add("Cluster/before/hTrackdEtadPhi", "d#eta vs. d#varphi of matched tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDEta, thAxisDPhi}, true);
  fRegistry->add("Cluster/before/hTrackdEtaPt", "d#eta vs. track pT of matched tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDEta, thAxisMomentum}, true);
  fRegistry->add("Cluster/before/hTrackdPhiPt", "d#varphi vs. track pT of matched tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDPhi, thAxisMomentum}, true);
  fRegistry->add("Cluster/before/hSecTrackdEtadPhi", "d#eta vs. d#varphi of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDEta, thAxisDPhi}, true);
  fRegistry->add("Cluster/before/hSecTrackdEtaPt", "d#eta vs. track pT of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDEta, thAxisMomentum}, true);
  fRegistry->add("Cluster/before/hSecTrackdPhiPt", "d#varphi vs. track pT of matched secondary tracks;d#eta;d#varphi (rad.)", o2::framework::kTH2F, {thAxisDPhi, thAxisMomentum}, true);

  if (do2DQA) { // Check if 2D QA histograms were selected in em-qc task
    fRegistry->add("Cluster/before/hNCell", "#it{N}_{cells};N_{cells} (GeV);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{26, -0.5, 25.5}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hM02", "Long ellipse axis;#it{M}_{02} (cm);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, 0, 2}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTime", "Cluster time;#it{t}_{cls} (ns);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{300, -150, 150}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdEta", "d#eta vs. E of matched tracks;d#eta;#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {thAxisDEta, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackdPhi", "d#phi vs. E of matched tracks;d#varphi (rad.);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {thAxisDPhi, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c});#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, 0., 5.}, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hSecTrackdEta", "d#eta vs. E of matched tracks;d#eta;#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {thAxisDEta, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hSecTrackdPhi", "d#phi vs. E of matched tracks;d#varphi (rad.);#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {thAxisDPhi, {200, 0, 20}}, true);
    fRegistry->add("Cluster/before/hSecTrackEOverP", "Energy of cluster divided by momentum of matched tracks;#it{E}_{cluster}/#it{p}_{track} (#it{c});#it{E}_{cluster} (GeV)", o2::framework::kTH2F, {{200, 0., 5.}, {200, 0, 20}}, true);
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

inline void fillTrackQA1D(
  o2::framework::HistogramRegistry* reg,
  const auto& histbase,
  float e, float p, float w)
{
  reg->fill(histbase + HIST("hTrackEOverP"), e / p, w);
}

inline void fillTrackQA2D(
  o2::framework::HistogramRegistry* reg,
  const auto& histbase,
  float e, float dEta, float dPhi, float p, float pt, float w)
{
  reg->fill(histbase + HIST("hTrackdEtadPhi"), dEta, dPhi, w);
  reg->fill(histbase + HIST("hTrackEOverP"), e / p, e, w);
  reg->fill(histbase + HIST("hTrackdEta"), dEta, e, w);
  reg->fill(histbase + HIST("hTrackdPhi"), dPhi, e, w);
  reg->fill(histbase + HIST("hTrackdEtaPt"), dEta, pt, w);
  reg->fill(histbase + HIST("hTrackdPhiPt"), dPhi, pt, w);
}

template <const int cls_id, o2::soa::is_iterator TCluster, is_optional_table TMatchedTracks = std::nullptr_t, is_optional_table TMatchedSecondaries = std::nullptr_t>
inline void fillClusterHistograms(o2::framework::HistogramRegistry* fRegistry, TCluster cluster, bool do2DQA, float weight = 1.f, TMatchedTracks const& primTracks = nullptr, TMatchedSecondaries const& secTracks = nullptr)
{
  const auto e = cluster.e();
  const auto eta = cluster.eta();
  const auto phi = cluster.phi();

  static constexpr std::string_view ClusterTypes[2] = {"Cluster/before/", "Cluster/after/"};

  if constexpr (HasTrackMatching<TCluster>) {
    for (size_t iTrack = 0; iTrack < cluster.deltaEta().size(); ++iTrack) {
      const float dEta = cluster.deltaEta()[iTrack];
      const float dPhi = cluster.deltaPhi()[iTrack];
      const float p = cluster.trackp()[iTrack];
      const float pt = cluster.trackpt()[iTrack];

      if (do2DQA) {
        fillTrackQA2D(fRegistry, HIST(ClusterTypes[cls_id]), e, dEta, dPhi, p, pt, weight);
      } else {
        fillTrackQA1D(fRegistry, HIST(ClusterTypes[cls_id]), e, p, weight);
      }
    }
  }

  if constexpr (HasPrimaries<TMatchedTracks> && IsTrackContainer<TMatchedTracks>) {
    for (const auto& primTrack : primTracks) {
      if (do2DQA) {
        fillTrackQA2D(fRegistry, HIST(ClusterTypes[cls_id]), e, primTrack.deltaEta(), primTrack.deltaPhi(), primTrack.trackP(), primTrack.trackPt(), weight);
      } else {
        fillTrackQA1D(fRegistry, HIST(ClusterTypes[cls_id]), e, primTrack.trackP(), weight);
      }
    }
  }

  if constexpr (HasSecondaries<TMatchedSecondaries> && IsTrackContainer<TMatchedSecondaries>) {
    for (const auto& secTrack : secTracks) {
      const auto trackP = secTrack.trackP();
      if (do2DQA) {
        const auto dEta = secTrack.deltaEta();
        const auto dPhi = secTrack.deltaPhi();
        const auto trackPt = secTrack.trackPt();
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackdEtadPhi"), dEta, dPhi, weight);
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackEOverP"), e / trackP, e, weight);
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackdEta"), dEta, e, weight);
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackdPhi"), dPhi, e, weight);
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackdEtaPt"), dEta, trackPt, weight);
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackdPhiPt"), dPhi, trackPt, weight);
      } else {
        fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hSecTrackEOverP"), e / trackP, weight);
      }
    }
  }

  fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hE"), e, weight);
  fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hPt"), cluster.pt(), weight);
  fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hEtaPhi"), eta, phi, weight);
  if constexpr (HasTrackMatching<TCluster>) {
    // primary matched tracks are stored directly inside cluster
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hNTracks"),
                    cluster.deltaEta().size(), weight);
  } else if constexpr (HasPrimaries<TMatchedTracks> && IsTrackContainer<TMatchedTracks>) {
    // primary matched tracks are stored as their own table
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hNTracks"), primTracks.size(), weight);
  }
  if constexpr (HasSecondaries<TMatchedSecondaries> && IsTrackContainer<TMatchedSecondaries>) {
    // secondary matched tracks are stored as their own table
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hNSecTracks"), secTracks.size(), weight);
  }

  if (do2DQA) {
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hNCell"), cluster.nCells(), e, weight);
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hM02"), cluster.m02(), e, weight);
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hTime"), cluster.time(), e, weight);
  } else {
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hNCell"), cluster.nCells(), weight);
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hM02"), cluster.m02(), weight);
    fRegistry->fill(HIST(ClusterTypes[cls_id]) + HIST("hTime"), cluster.time(), weight);
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::clusterhistogram

#endif // PWGEM_PHOTONMESON_UTILS_CLUSTERHISTOGRAMS_H_
