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

/// \file EmcalMatchingUtilities.h
/// \brief EMCal track matching utilities
///
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>, Goethe-University Frankfurt

#ifndef PWGJE_CORE_EMCALMATCHINGUTILITIES_H_
#define PWGJE_CORE_EMCALMATCHINGUTILITIES_H_

#include "Common/Core/RecoDecay.h"

#include "CommonConstants/MathConstants.h"

#include <TKDTree.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/parameters.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility> // for std::pair
#include <vector>

namespace emcmatchingutilities
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

/**
 * Match cells to tracks.
 *
 * Match cells to tracks, where maxNumberMatches are considered in dR=maxMatchingDistance.
 * If no unique cell match was found for a track, an index of -1 is stored..
 *
 * @param cellPhi cell collection phi.
 * @param cellEta cell collection eta.
 * @param trackPhi track collection phi.
 * @param trackEta track collection eta.
 * @param maxMatchingDistance Maximum matching distance.
 *
 * @returns (track to cell index map)
 */
template <typename T>
void matchCellsAndTracks(
  std::vector<T>& cellPhi,
  std::vector<T>& cellEta,
  std::vector<T>& trackPhi,
  std::vector<T>& trackEta,
  double maxMatchingDistance,
  std::vector<int>& matchIndexCell)
{
  // Validation
  const std::size_t nClusters = cellEta.size();
  const std::size_t nTracks = trackEta.size();

  if (!(nClusters && nTracks)) {
    // There are no jets, so nothing to be done.
    return;
  }
  // Input sizes must match
  if (cellPhi.size() != cellEta.size()) {
    throw std::invalid_argument("Cells collection eta and phi sizes don't match! Check the inputs!");
    return;
  }
  if (trackPhi.size() != trackEta.size()) {
    throw std::invalid_argument("Track collection eta and phi sizes don't match! Check the inputs!");
    return;
  }

  matchIndexCell.assign(nTracks, -1);

  // Build the KD-trees using vectors
  // We build two trees:
  // treeBase, which contains the base collection.
  // treeTag, which contains the tag collection.
  // The trees are built to match in two dimensions (eta, phi)
  TKDTree<int, T> treeCell(cellEta.size(), 2, 1);
  // By utilizing SetData, we can avoid having to copy the data again.
  treeCell.SetData(0, cellEta.data());
  treeCell.SetData(1, cellPhi.data());
  treeCell.Build();

  // Find the cell closest to each track
  for (std::size_t iTrack = 0; iTrack < nTracks; iTrack++) {
    T point[2] = {trackEta[iTrack], trackPhi[iTrack]};
    int index[1];  // size 50 for safety
    T distance[1]; // size 50 for safery
    std::fill_n(index, 1, -1);
    std::fill_n(distance, 1, std::numeric_limits<T>::max());
    treeCell.FindNearestNeighbors(point, 1, index, distance);
    if (index[0] >= 0 && distance[0] < maxMatchingDistance) {
      matchIndexCell[iTrack] = index[0];
    }
  }
  return;
}

template <typename T>
void matchCellsAndTracks2(
  std::vector<T>& cellPhi,
  std::vector<T>& cellEta,
  std::vector<T>& trackPhi,
  std::vector<T>& trackEta,
  double maxMatchingDistance,
  std::vector<int>& matchIndexCell)
{
  using Point = boost::geometry::model::point<T, 2, boost::geometry::cs::cartesian>;
  using Value = std::pair<Point, int>;

  const std::size_t nCells = cellEta.size();
  const std::size_t nTracks = trackEta.size();

  if (!(nCells && nTracks)) {
    return;
  }

  if (cellPhi.size() != nCells) {
    throw std::invalid_argument("Cells collection eta and phi sizes don't match! Check the inputs!");
  }
  if (trackPhi.size() != nTracks) {
    throw std::invalid_argument("Track collection eta and phi sizes don't match! Check the inputs!");
  }

  // Build R-tree from cells
  std::vector<Value> rtreeEntries;
  rtreeEntries.reserve(nCells);

  for (std::size_t i = 0; i < nCells; ++i) {
    rtreeEntries.emplace_back(Point(cellEta[i], cellPhi[i]), static_cast<int>(i));
  }

  boost::geometry::index::rtree<Value, boost::geometry::index::quadratic<16>> rtree(rtreeEntries);

  matchIndexCell.assign(nTracks, -1);

  for (std::size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    Point query(trackEta[iTrack], trackPhi[iTrack]);

    std::vector<Value> result;
    rtree.query(boost::geometry::index::nearest(query, 1), std::back_inserter(result));

    if (!result.empty()) {
      const auto& [matchedPoint, matchedIdx] = result.front();

      // Compute actual distance squared
      T dEta = trackEta[iTrack] - boost::geometry::get<0>(matchedPoint);
      T dPhiRaw = trackPhi[iTrack] - boost::geometry::get<1>(matchedPoint);
      T dPhi = std::fabs(dPhiRaw);
      if (dPhi > static_cast<T>(M_PI)) {
        dPhi = static_cast<T>(2 * M_PI) - dPhi;
      }
      T dist2 = dEta * dEta + dPhi * dPhi;

      if (dist2 < maxMatchingDistance * maxMatchingDistance) {
        matchIndexCell[iTrack] = matchedIdx;
      }
    }
  }
}

template <typename T, typename U>
float deltaR(T const& A, U const& B)
{
  float dPhi = RecoDecay::constrainAngle(A.phi() - B.phi(), -o2::constants::math::PI);
  float dEta = A.eta() - B.eta();

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}
// same as deltaR but explicit specification of the eta and phi components
template <typename T, typename U, typename V, typename W>
float deltaR(T const& eta1, U const& phi1, V const& eta2, W const& phi2)
{
  float dPhi = RecoDecay::constrainAngle(phi1 - phi2, -o2::constants::math::PI);
  float dEta = eta1 - eta2;

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}
}; // namespace emcmatchingutilities

#endif // PWGJE_CORE_EMCALMATCHINGUTILITIES_H_
