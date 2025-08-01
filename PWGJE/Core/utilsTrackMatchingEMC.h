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

/// \file utilsTrackMatchingEMC.h
/// \brief EMCal track matching related utils
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>

#ifndef PWGJE_CORE_UTILSTRACKMATCHINGEMC_H_
#define PWGJE_CORE_UTILSTRACKMATCHINGEMC_H_

#include <TKDTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace tmemcutilities
{

struct MatchResult {
  std::vector<std::vector<int>> matchIndexTrack;
  std::vector<std::vector<float>> matchDeltaPhi;
  std::vector<std::vector<float>> matchDeltaEta;
};

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
MatchResult matchTracksToCluster(
  std::span<float> clusterPhi,
  std::span<float> clusterEta,
  std::span<float> trackPhi,
  std::span<float> trackEta,
  double maxMatchingDistance,
  int maxNumberMatches)
{
  const std::size_t nClusters = clusterEta.size();
  const std::size_t nTracks = trackEta.size();
  MatchResult result;

  if (nClusters == 0 || nTracks == 0) {
    // There are no jets, so nothing to be done.
    return result;
  }
  // Input sizes must match
  if (clusterPhi.size() != clusterEta.size()) {
    throw std::invalid_argument("cluster collection eta and phi sizes don't match. Check the inputs.");
  }
  if (trackPhi.size() != trackEta.size()) {
    throw std::invalid_argument("track collection eta and phi sizes don't match. Check the inputs.");
  }

  result.matchIndexTrack.resize(nClusters);
  result.matchDeltaPhi.resize(nClusters);
  result.matchDeltaEta.resize(nClusters);

  // Build the KD-trees using vectors
  // We build two trees:
  // treeBase, which contains the base collection.
  // treeTag, which contains the tag collection.
  // The trees are built to match in two dimensions (eta, phi)
  TKDTree<int, float> treeTrack(trackEta.size(), 2, 1);
  treeTrack.SetData(0, trackEta.data());
  treeTrack.SetData(1, trackPhi.data());
  treeTrack.Build();

  // Find the track closest to each cluster.
  for (std::size_t iCluster = 0; iCluster < nClusters; iCluster++) {
    float point[2] = {clusterEta[iCluster], clusterPhi[iCluster]};
    int index[50];      // size 50 for safety
    float distance[50]; // size 50 for safery
    std::fill_n(index, 50, -1);
    std::fill_n(distance, 50, std::numeric_limits<float>::max());
    treeTrack.FindNearestNeighbors(point, maxNumberMatches, index, distance);

    // allocate enough memory
    result.matchIndexTrack[iCluster].reserve(maxNumberMatches);
    result.matchDeltaPhi[iCluster].reserve(maxNumberMatches);
    result.matchDeltaEta[iCluster].reserve(maxNumberMatches);

    // test whether indices are matching:
    for (int m = 0; m < maxNumberMatches; m++) {
      if (index[m] >= 0 && distance[m] < maxMatchingDistance) {
        result.matchIndexTrack[iCluster].push_back(index[m]);
        result.matchDeltaPhi[iCluster].push_back(trackPhi[index[m]] - clusterPhi[iCluster]);
        result.matchDeltaEta[iCluster].push_back(trackEta[index[m]] - clusterEta[iCluster]);
      }
    }
  }
  return result;
}
}; // namespace tmemcutilities

#endif // PWGJE_CORE_UTILSTRACKMATCHINGEMC_H_
