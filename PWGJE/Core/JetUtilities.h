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

#ifndef PWGJE_CORE_JETUTILITIES_H_
#define PWGJE_CORE_JETUTILITIES_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

#include <TKDTree.h>

#include "Framework/Logger.h"

namespace JetUtilities
{
/**
 * Duplicates jets around the phi boundary which are within the matching distance.
 *
 * NOTE: Assumes, but does not validate, that 0 <= phi < 2pi.
 *
 * @param jetsPhi Jets phi
 * @param jetsEta Jets eta
 * @param maxMatchingDistance Maximum matching distance. Only duplicate jets within this distance of the boundary.
 */
template <typename T>
std::tuple<std::vector<std::size_t>, std::vector<T>, std::vector<T>> DuplicateJetsAroundPhiBoundary(
  std::vector<T>& jetsPhi,
  std::vector<T>& jetsEta,
  double maxMatchingDistance,
  // TODO: Remove additional margin after additional testing.
  double additionalMargin = 0.05)
{
  const std::size_t nJets = jetsPhi.size();
  std::vector<std::size_t> jetsMapToJetIndex(nJets);
  // We need to keep track of the map from the duplicated vector to the existing jets.
  // To start, we fill this map (practically, it maps from vector index to an index, so we
  // just use a standard vector) with the existing jet indices, which range from 0..nJets.
  std::iota(jetsMapToJetIndex.begin(), jetsMapToJetIndex.end(), 0);

  // The full duplicated jets will be stored in a new vector to avoid disturbing the input data.
  std::vector<T> jetsPhiComparison(jetsPhi);
  std::vector<T> jetsEtaComparison(jetsEta);

  // Duplicate the jets
  // When a jet is outside of the desired phi range, we make a copy of the jet, shifting it by 2pi
  // as necessary. The eta value is copied directly, as is the index of the original jet, so we can
  // map from index in duplicated data -> index in original collection.
  // NOTE: Assumes, but does not validate, that 0 <= phi < 2pi.
  for (std::size_t i = 0; i < nJets; i++) {
    // Handle lower edge
    if (jetsPhi[i] <= (maxMatchingDistance + additionalMargin)) {
      jetsPhiComparison.emplace_back(jetsPhi[i] + 2 * M_PI);
      jetsEtaComparison.emplace_back(jetsEta[i]);
      jetsMapToJetIndex.emplace_back(jetsMapToJetIndex[i]);
    }
    // Handle upper edge
    if (jetsPhi[i] >= (2 * M_PI - (maxMatchingDistance + additionalMargin))) {
      jetsPhiComparison.emplace_back(jetsPhi[i] - 2 * M_PI);
      jetsEtaComparison.emplace_back(jetsEta[i]);
      jetsMapToJetIndex.emplace_back(jetsMapToJetIndex[i]);
    }
  }

  return {jetsMapToJetIndex, jetsPhiComparison, jetsEtaComparison};
}

/**
 * Implementation of geometrical jet matching.
 *
 * Jets are required to match uniquely - namely: base <-> tag. Only one direction of matching isn't enough.
 * Unless special conditions are required, it's better to use `MatchJetsGeometrically`, which has an
 * easier to use interface.
 *
 * NOTE: The vectors for matching could all be const, except that SetData() doesn't take a const.
 *
 * @param jetsBasePhi Base jet collection phi.
 * @param jetsBaseEta Base jet collection eta.
 * @param jetsBasePhiForMatching Base jet collection phi to use for matching.
 * @param jetsBaseEtaForMatching Base jet collection eta to use for matching.
 * @param jetMapBaseToJetIndex Base jet collection index map from duplicated jets to original jets.
 * @param jetsTagPhi Tag jet collection phi.
 * @param jetsTagEta Tag jet collection eta.
 * @param jetsTagPhiForMatching Tag jet collection phi to use for matching.
 * @param jetsTagEtaForMatching Tag jet collection eta to use for matching.
 * @param jetMapTagToJetIndex Tag jet collection index map from duplicated jets to original jets.
 * @param maxMatchingDistance Maximum matching distance.
 *
 * @returns (Base to tag index map, tag to base index map) for uniquely matched jets.
 */
template <typename T>
std::tuple<std::vector<int>, std::vector<int>> MatchJetsGeometricallyImpl(
  const std::vector<T>& jetsBasePhi,
  const std::vector<T>& jetsBaseEta,
  std::vector<T> jetsBasePhiForMatching,
  std::vector<T> jetsBaseEtaForMatching,
  const std::vector<std::size_t> jetMapBaseToJetIndex,
  const std::vector<T>& jetsTagPhi,
  const std::vector<T>& jetsTagEta,
  std::vector<T> jetsTagPhiForMatching,
  std::vector<T> jetsTagEtaForMatching,
  const std::vector<std::size_t> jetMapTagToJetIndex,
  double maxMatchingDistance)
{
  // Validation
  // If no jets in either collection, then return immediately.
  const std::size_t nJetsBase = jetsBaseEta.size();
  const std::size_t nJetsTag = jetsTagEta.size();
  if (!(nJetsBase && nJetsTag)) {
    return std::make_tuple(std::vector<int>(nJetsBase, -1), std::vector<int>(nJetsTag, -1));
  }
  // Require that the comparison vectors are greater than or equal to the standard collections.
  if (jetsBasePhiForMatching.size() < jetsBasePhi.size()) {
    throw std::invalid_argument("Base collection phi for matching is smaller than the input base collection.");
  }
  if (jetsBaseEtaForMatching.size() < jetsBaseEta.size()) {
    throw std::invalid_argument("Base collection eta for matching is smaller than the input base collection.");
  }
  if (jetsTagPhiForMatching.size() < jetsTagPhi.size()) {
    throw std::invalid_argument("Tag collection phi for matching is smaller than the input tag collection.");
  }
  if (jetsTagEtaForMatching.size() < jetsTagEta.size()) {
    throw std::invalid_argument("Tag collection eta for matching is smaller than the input tag collection.");
  }

  // Build the KD-trees using vectors
  // We build two trees:
  // treeBase, which contains the base collection.
  // treeTag, which contains the tag collection.
  // The trees are built to match in two dimensions (eta, phi)
  TKDTree<int, T> treeBase(jetsBaseEtaForMatching.size(), 2, 1), treeTag(jetsTagEtaForMatching.size(), 2, 1);
  // By utilizing SetData, we can avoid having to copy the data again.
  treeBase.SetData(0, jetsBaseEtaForMatching.data());
  treeBase.SetData(1, jetsBasePhiForMatching.data());
  treeBase.Build();
  treeTag.SetData(0, jetsTagEtaForMatching.data());
  treeTag.SetData(1, jetsTagPhiForMatching.data());
  treeTag.Build();

  // Storage for the jet matching indices.
  // matchIndexTag maps from the base index to the tag index.
  // matchBaseTag maps from the tag index to the base index.
  std::vector<int> matchIndexTag(nJetsBase, -1), matchIndexBase(nJetsTag, -1);

  // Find the tag jet closest to each base jet.
  for (std::size_t iBase = 0; iBase < nJetsBase; iBase++) {
    Double_t point[2] = {jetsBaseEta[iBase], jetsBasePhi[iBase]};
    Int_t index(-1);
    Double_t distance(-1);
    treeTag.FindNearestNeighbors(point, 1, &index, &distance);
    // test whether indices are matching:
    if (index >= 0 && distance < maxMatchingDistance) {
      LOG(debug) << "Found closest tag jet for " << iBase << " with match index " << index << " and distance " << distance << "\n";
      matchIndexTag[iBase] = index;
    } else {
      LOG(debug) << "Closest tag jet not found for " << iBase << ", distance to closest " << distance << "\n";
    }
  }

  // Find the base jet closest to each tag jet
  for (std::size_t iTag = 0; iTag < nJetsTag; iTag++) {
    Double_t point[2] = {jetsTagEta[iTag], jetsTagPhi[iTag]};
    Int_t index(-1);
    Double_t distance(-1);
    treeBase.FindNearestNeighbors(point, 1, &index, &distance);
    if (index >= 0 && distance < maxMatchingDistance) {
      LOG(debug) << "Found closest base jet for " << iTag << " with match index " << index << " and distance " << distance << std::endl;
      matchIndexBase[iTag] = index;
    } else {
      LOG(debug) << "Closest tag jet not found for " << iTag << ", distance to closest " << distance << "\n";
    }
  }

  // Convert indices in the duplicated jet vectors into the original jet indices.
  // First for the base -> tag map.
  for (auto& v : matchIndexTag) {
    // If it's -1, it means that it didn't find a matching jet.
    // We have to explicitly check for it here because it would be an invalid index.
    if (v != -1) {
      v = jetMapTagToJetIndex[v];
    }
  }
  // Then for the index -> base map.
  for (auto& v : matchIndexBase) {
    if (v != -1) {
      v = jetMapBaseToJetIndex[v];
    }
  }

  // Finally, we'll check for true matches, which are pairs where the base jet is the
  // closest to the tag jet and vice versa
  // As the lists are linear a loop over the outer base jet is sufficient.
  std::vector<int> baseToTagMap(nJetsBase, -1);
  std::vector<int> tagToBaseMap(nJetsTag, -1);
  LOG(debug) << "Starting true jet loop: nbase(" << nJetsBase << "), ntag(" << nJetsTag << ")\n";
  for (std::size_t iBase = 0; iBase < nJetsBase; iBase++) {
    LOG(debug) << "base jet " << iBase << ": match index in tag jet container " << matchIndexTag[iBase] << "\n";
    if (matchIndexTag[iBase] > -1) {
      LOG(debug) << "tag jet " << matchIndexTag[iBase] << ": matched base jet " << matchIndexBase[matchIndexTag[iBase]] << "\n";
    }
    if (matchIndexTag[iBase] > -1 && matchIndexBase[matchIndexTag[iBase]] == static_cast<int>(iBase)) {
      LOG(debug) << "True match! base index: " << iBase << ", tag index: " << matchIndexTag[iBase] << "\n";
      baseToTagMap[iBase] = matchIndexTag[iBase];
      tagToBaseMap[matchIndexTag[iBase]] = iBase;
    }
  }

  return std::make_tuple(baseToTagMap, tagToBaseMap);
}

/**
 * Geometrical jet matching.
 *
 * Match jets in the "base" collection with those in the "tag" collection. Jets are matched within
 * the provided matching distance. Jets are required to match uniquely - namely: base <-> tag.
 * Only one direction of matching isn't enough.
 *
 * If no unique match was found for a jet, an index of -1 is stored.
 *
 * @param jetsBasePhi Base jet collection phi.
 * @param jetsBaseEta Base jet collection eta.
 * @param jetsTagPhi Tag jet collection phi.
 * @param jetsTagEta Tag jet collection eta.
 * @param maxMatchingDistance Maximum matching distance.
 *
 * @returns (Base to tag index map, tag to base index map) for uniquely matched jets.
 */
template <typename T>
std::tuple<std::vector<int>, std::vector<int>> MatchJetsGeometrically(
  std::vector<T> jetsBasePhi,
  std::vector<T> jetsBaseEta,
  std::vector<T> jetsTagPhi,
  std::vector<T> jetsTagEta,
  double maxMatchingDistance)
{
  // Validation
  const std::size_t nJetsBase = jetsBaseEta.size();
  const std::size_t nJetsTag = jetsTagEta.size();
  if (!(nJetsBase && nJetsTag)) {
    // There are no jets, so nothing to be done.
    return std::make_tuple(std::vector<int>(nJetsBase, -1), std::vector<int>(nJetsTag, -1));
  }
  // Input sizes must match
  if (jetsBasePhi.size() != jetsBaseEta.size()) {
    throw std::invalid_argument("Base collection eta and phi sizes don't match. Check the inputs.");
  }
  if (jetsTagPhi.size() != jetsTagEta.size()) {
    throw std::invalid_argument("Tag collection eta and phi sizes don't match. Check the inputs.");
  }

  // To perform matching with periodic boundary conditions (ie. phi) with a KDTree, we need
  // to duplicate data up to maxMatchingDistance in phi because phi is periodic.
  // NOTE: vectors are modified in place to avoid copies.
  auto&& [jetMapBaseToJetIndex, jetsBasePhiComparison, jetsBaseEtaComparison] = DuplicateJetsAroundPhiBoundary(jetsBasePhi, jetsBaseEta, maxMatchingDistance);
  auto&& [jetMapTagToJetIndex, jetsTagPhiComparison, jetsTagEtaComparison] = DuplicateJetsAroundPhiBoundary(jetsTagPhi, jetsTagEta, maxMatchingDistance);

  // Finally, perform the actual matching.
  auto&& [baseToTagMap, tagToBaseMap] = MatchJetsGeometricallyImpl(
    jetsBasePhi, jetsBaseEta, jetsBasePhiComparison, jetsBaseEtaComparison, jetMapBaseToJetIndex,
    jetsTagPhi, jetsTagEta, jetsTagPhiComparison, jetsTagEtaComparison, jetMapTagToJetIndex,
    maxMatchingDistance);

  return std::make_tuple(baseToTagMap, tagToBaseMap);
}

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

  // for (std::size_t iTrack = 0; iTrack < nTracks; iTrack++) {
  //   if (trackEta[iTrack] == 0)
  //     LOG(warning) << "Track eta is 0!";
  // }

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

}; // namespace JetUtilities

#endif // PWGJE_CORE_JETUTILITIES_H_
