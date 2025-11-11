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

/// \file JetUtilities.h
/// \brief Jet related utilities
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETUTILITIES_H_
#define PWGJE_CORE_JETUTILITIES_H_

#include "Common/Core/RecoDecay.h"

#include <TKDTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <math.h>

namespace jetutilities
{

/**
 * Match clusters and tracks.
 *
 * Match cluster with tracks, where maxNumberMatches are considered in dR=maxMatchingDistance.
 * If no unique match was found for a jet, an index of -1 is stored.
 * The same map is created for clusters matched to tracks e.g. for electron analyses.
 *
 * @param clusterPhi cluster collection phi.
 * @param clusterEta cluster collection eta.
 * @param trackPhi track collection phi.
 * @param trackEta track collection eta.
 * @param maxMatchingDistance Maximum matching distance.
 * @param maxNumberMatches Maximum number of matches (e.g. 5 closest).
 *
 * @returns (cluster to track index map, track to cluster index map)
 */
template <typename T>
std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> MatchClustersAndTracks(
  std::vector<T>& clusterPhi,
  std::vector<T>& clusterEta,
  std::vector<T>& trackPhi,
  std::vector<T>& trackEta,
  double maxMatchingDistance,
  int maxNumberMatches)
{
  // test
  // Validation
  const std::size_t nClusters = clusterEta.size();
  const std::size_t nTracks = trackEta.size();
  if (!(nClusters && nTracks)) {
    // There are no jets, so nothing to be done.
    return std::make_tuple(std::vector<std::vector<int>>(nClusters, std::vector<int>(maxNumberMatches, -1)), std::vector<std::vector<int>>(nTracks, std::vector<int>(maxNumberMatches, -1)));
  }
  // Input sizes must match
  if (clusterPhi.size() != clusterEta.size()) {
    throw std::invalid_argument("cluster collection eta and phi sizes don't match. Check the inputs.");
  }
  if (trackPhi.size() != trackEta.size()) {
    throw std::invalid_argument("track collection eta and phi sizes don't match. Check the inputs.");
  }

  // Build the KD-trees using vectors
  // We build two trees:
  // treeBase, which contains the base collection.
  // treeTag, which contains the tag collection.
  // The trees are built to match in two dimensions (eta, phi)
  TKDTree<int, T> treeCluster(clusterEta.size(), 2, 1), treeTrack(trackEta.size(), 2, 1);
  // By utilizing SetData, we can avoid having to copy the data again.
  treeCluster.SetData(0, clusterEta.data());
  treeCluster.SetData(1, clusterPhi.data());
  treeCluster.Build();
  treeTrack.SetData(0, trackEta.data());
  treeTrack.SetData(1, trackPhi.data());
  treeTrack.Build();

  // Storage for the cluster matching indices.
  std::vector<std::vector<int>> matchIndexTrack(nClusters, std::vector<int>(maxNumberMatches, -1));
  std::vector<std::vector<int>> matchIndexCluster(nTracks, std::vector<int>(maxNumberMatches, -1));

  // Find the track closest to each cluster.
  for (std::size_t iCluster = 0; iCluster < nClusters; iCluster++) {
    T point[2] = {clusterEta[iCluster], clusterPhi[iCluster]};
    int index[50];  // size 50 for safety
    T distance[50]; // size 50 for safery
    std::fill_n(index, 50, -1);
    std::fill_n(distance, 50, std::numeric_limits<T>::max());
    treeTrack.FindNearestNeighbors(point, maxNumberMatches, index, distance);
    // test whether indices are matching:
    matchIndexTrack[iCluster] = std::vector<int>(maxNumberMatches);
    for (int m = 0; m < maxNumberMatches; m++) {
      if (index[m] >= 0 && distance[m] < maxMatchingDistance) {
        matchIndexTrack[iCluster][m] = index[m];
      } else {
        // no match or no more matches found, fill -1
        matchIndexTrack[iCluster][m] = -1;
      }
    }
  }

  // Find the base jet closest to each tag jet
  for (std::size_t iTrack = 0; iTrack < nTracks; iTrack++) {
    T point[2] = {trackEta[iTrack], trackPhi[iTrack]};
    int index[50];  // size 50 for safety
    T distance[50]; // size 50 for safery
    std::fill_n(index, 50, -1);
    std::fill_n(distance, 50, std::numeric_limits<T>::max());
    treeCluster.FindNearestNeighbors(point, maxNumberMatches, index, distance);
    matchIndexCluster[iTrack] = std::vector<int>(maxNumberMatches);
    // loop over maxNumberMatches closest matches
    for (int m = 0; m < maxNumberMatches; m++) {
      if (index[m] >= 0 && distance[m] < maxMatchingDistance) {
        matchIndexCluster[iTrack][m] = index[m];
      } else {
        // no match jet or no more matches found, fill -1
        matchIndexCluster[iTrack][m] = -1;
      }
    }
  }
  return std::make_tuple(matchIndexTrack, matchIndexCluster);
}

template <typename T, typename U>
float deltaR(T const& A, U const& B)
{
  float dPhi = RecoDecay::constrainAngle(A.phi() - B.phi(), -M_PI);
  float dEta = A.eta() - B.eta();

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}
// same as deltaR but explicit specification of the eta and phi components
template <typename T, typename U, typename V, typename W>
float deltaR(T const& eta1, U const& phi1, V const& eta2, W const& phi2)
{
  float dPhi = RecoDecay::constrainAngle(phi1 - phi2, -M_PI);
  float dEta = eta1 - eta2;

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

/// @brief Background estimator using the perpendicular cone method
/// @param inputParticles
/// @param jet
/// @return Rho, RhoM the underlying event density

template <typename T, typename U, typename V>
std::tuple<double, double> estimateRhoPerpCone(const T& inputParticles, const U& jet, V perpConeR)
{

  if (inputParticles.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  double perpPtDensity1 = 0;
  double perpPtDensity2 = 0;
  double perpMdDensity1 = 0;
  double perpMdDensity2 = 0;

  const double jetPhi = RecoDecay::constrainAngle<double, double>(jet.phi(), -M_PI);
  const double jetEta = jet.eta();
  const double radius = static_cast<double>(perpConeR);

  // build 2 perp cones in phi around the leading jet (right and left of the jet)
  double PerpendicularConeAxisPhi1 = RecoDecay::constrainAngle<double, double>(jetPhi + (M_PI / 2.), -M_PI); // This will contrain the angel between -pi & Pi
  double PerpendicularConeAxisPhi2 = RecoDecay::constrainAngle<double, double>(jetPhi - (M_PI / 2.), -M_PI); // This will contrain the angel between -pi & Pi

  for (const auto& particle : inputParticles) {
    // sum the momentum of all paricles that fill the two cones
    const double phi = RecoDecay::constrainAngle<double, double>(particle.phi(), -M_PI);
    double dPhi1 = RecoDecay::constrainAngle<double, double>(phi - PerpendicularConeAxisPhi1, -M_PI); // This will contrain the angel between -pi & Pi
    double dPhi2 = RecoDecay::constrainAngle<double, double>(phi - PerpendicularConeAxisPhi2, -M_PI); // This will contrain the angel between -pi & Pi
    double dEta = jetEta - particle.eta();                                                            // The perp cone eta is the same as the leading jet since the cones are perpendicular only in phi
    if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= static_cast<double>(radius)) {
      perpPtDensity1 += particle.pt();
      perpMdDensity1 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }

    if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= static_cast<double>(radius)) {
      perpPtDensity2 += particle.pt();
      perpMdDensity2 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }
  }

  // Caculate rho as the ratio of average pT of the two cones / the cone area
  double perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * M_PI * static_cast<double>(radius) * static_cast<double>(radius));
  double perpMdDensity = (perpMdDensity1 + perpMdDensity2) / (2 * M_PI * static_cast<double>(radius) * static_cast<double>(radius));

  return std::make_tuple(perpPtDensity, perpMdDensity);
}

}; // namespace jetutilities

#endif // PWGJE_CORE_JETUTILITIES_H_
