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

/// \file JetMatchingUtilities.h
/// \brief Jet Matching related utilities
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Aimeric Lanodu <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETMATCHINGUTILITIES_H_
#define PWGJE_CORE_JETMATCHINGUTILITIES_H_

#include "PWGJE/Core/JetCandidateUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/Logger.h>

#include <TKDTree.h>

#include <RtypesCore.h>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <math.h>

namespace jetmatchingutilities
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
  const std::vector<std::size_t>& jetMapBaseToJetIndex,
  const std::vector<T>& jetsTagPhi,
  const std::vector<T>& jetsTagEta,
  std::vector<T> jetsTagPhiForMatching,
  std::vector<T> jetsTagEtaForMatching,
  const std::vector<std::size_t>& jetMapTagToJetIndex,
  const double maxMatchingDistance)
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

template <typename T, typename U>
void MatchGeo(T const& jetsBasePerCollision, U const& jetsTagPerCollision, std::vector<std::vector<int>>& baseToTagMatchingGeo, std::vector<std::vector<int>>& tagToBaseMatchingGeo, float maxMatchingDistance)
{
  std::vector<double> jetsR;
  for (const auto& jetBase : jetsBasePerCollision) {
    if (std::find(jetsR.begin(), jetsR.end(), std::round(jetBase.r())) == jetsR.end()) {
      jetsR.push_back(std::round(jetBase.r()));
    }
  }
  for (const auto& jetTag : jetsTagPerCollision) {
    if (std::find(jetsR.begin(), jetsR.end(), std::round(jetTag.r())) == jetsR.end()) {
      jetsR.push_back(std::round(jetTag.r()));
    }
  }
  for (auto jetR : jetsR) {
    std::vector<double> jetsBasePhi;
    std::vector<double> jetsBaseEta;
    std::vector<int> baseToTagMatchingGeoIndex;
    std::vector<int> tagToBaseMatchingGeoIndex;

    for (const auto& jetBase : jetsBasePerCollision) {
      if (std::round(jetBase.r()) != std::round(jetR)) {
        continue;
      }
      jetsBasePhi.emplace_back(jetBase.phi());
      jetsBaseEta.emplace_back(jetBase.eta());
    }
    std::vector<double> jetsTagPhi;
    std::vector<double> jetsTagEta;
    for (const auto& jetTag : jetsTagPerCollision) {
      if (std::round(jetTag.r()) != std::round(jetR)) {
        continue;
      }
      jetsTagPhi.emplace_back(jetTag.phi());
      jetsTagEta.emplace_back(jetTag.eta());
    }
    std::tie(baseToTagMatchingGeoIndex, tagToBaseMatchingGeoIndex) = MatchJetsGeometrically(jetsBasePhi, jetsBaseEta, jetsTagPhi, jetsTagEta, maxMatchingDistance); // change max distnace to a function call
    int jetBaseIndex = 0;
    int jetTagIndex = 0;
    for (const auto& jetBase : jetsBasePerCollision) {
      if (std::round(jetBase.r()) != std::round(jetR)) {
        continue;
      }
      jetTagIndex = baseToTagMatchingGeoIndex[jetBaseIndex];
      if (jetTagIndex > -1 && jetTagIndex < jetsTagPerCollision.size()) {
        int jetTagGlobalIndex = jetsTagPerCollision.iteratorAt(jetTagIndex).globalIndex();
        baseToTagMatchingGeo[jetBase.globalIndex()].push_back(jetTagGlobalIndex);
      }
      jetBaseIndex++;
    }
    for (const auto& jetTag : jetsTagPerCollision) {
      if (std::round(jetTag.r()) != std::round(jetR)) {
        continue;
      }
      jetBaseIndex = tagToBaseMatchingGeoIndex[jetTagIndex];
      if (jetBaseIndex > -1 && jetBaseIndex < jetsBasePerCollision.size()) {
        int jetBaseGlobalIndex = jetsBasePerCollision.iteratorAt(jetBaseIndex).globalIndex();
        tagToBaseMatchingGeo[jetTag.globalIndex()].push_back(jetBaseGlobalIndex);
      }
      jetTagIndex++;
    }
  }
}

// function that does the HF matching of jets from jetsBasePerColl and jets from jetsTagPerColl; assumes both jetsBasePerColl and jetsTagPerColl have access to Mc information
template <bool jetsBaseIsMc, bool jetsTagIsMc, typename T, typename U, typename V, typename M, typename N, typename O>
void MatchHF(T const& jetsBasePerCollision, U const& jetsTagPerCollision, std::vector<std::vector<int>>& baseToTagMatchingHF, std::vector<std::vector<int>>& tagToBaseMatchingHF, V const& /*candidatesBase*/, M const& /*candidatesTag*/, N const& tracksBase, O const& tracksTag)
{
  for (const auto& jetBase : jetsBasePerCollision) {
    if (jetBase.candidatesIds().size() == 0) {
      continue;
    }
    const auto candidateBase = jetBase.template candidates_first_as<V>();
    for (const auto& jetTag : jetsTagPerCollision) {
      if (jetTag.candidatesIds().size() == 0) {
        continue;
      }
      if (std::round(jetBase.r()) != std::round(jetTag.r())) {
        continue;
      }
      if constexpr (jetsBaseIsMc || jetsTagIsMc) {
        if (jetcandidateutilities::isMatchedCandidate(candidateBase)) {
          const auto candidateBaseMcId = jetcandidateutilities::matchedParticleId(candidateBase, tracksBase, tracksTag);
          const auto candidateTag = jetTag.template candidates_first_as<M>();
          const auto candidateTagId = candidateTag.mcParticleId();
          if (candidateBaseMcId == candidateTagId) {
            baseToTagMatchingHF[jetBase.globalIndex()].push_back(jetTag.globalIndex());
            tagToBaseMatchingHF[jetTag.globalIndex()].push_back(jetBase.globalIndex());
          }
        }
      } else {
        const auto candidateTag = jetTag.template candidates_first_as<M>();
        if (candidateBase.globalIndex() == candidateTag.globalIndex()) {
          baseToTagMatchingHF[jetBase.globalIndex()].push_back(jetTag.globalIndex());
          tagToBaseMatchingHF[jetTag.globalIndex()].push_back(jetBase.globalIndex());
        }
      }
    }
  }
}

template <bool isMc, typename T>
auto constexpr getConstituentId(T const& track)
{
  if constexpr (isMc) {

    if (track.has_mcParticle()) {
      return track.mcParticleId();
    } else {
      return -1;
    }

  } else {
    return track.globalIndex();
  }
}

template <bool isEMCAL, bool isCandidate, bool jetsBaseIsMc, bool jetsTagIsMc, typename T, typename U, typename V, typename O, typename P, typename Q, typename R, typename S>
float getPtSum(T const& tracksBase, U const& candidatesBase, V const& clustersBase, O const& tracksTag, P const& candidatesTag, Q const& clustersTag, R const& fullTracksBase, S const& fullTracksTag)
{
  std::vector<int> particleTracker;
  float ptSum = 0.;
  for (const auto& trackBase : tracksBase) {
    auto trackBaseId = getConstituentId<jetsTagIsMc>(trackBase);
    for (const auto& trackTag : tracksTag) {
      auto trackTagId = getConstituentId<jetsBaseIsMc>(trackTag);
      if (trackBaseId != -1 && trackBaseId == trackTagId) {
        ptSum += trackBase.pt();
        if constexpr (jetsBaseIsMc) {
          particleTracker.push_back(trackBaseId);
        }
        break;
      }
    }
  }
  if constexpr (isEMCAL) {
    if constexpr (jetsTagIsMc) {
      for (const auto& clusterBase : clustersBase) {
        for (const auto& clusterBaseParticleId : clusterBase.mcParticlesIds()) {
          bool isClusterMatched = false;
          for (const auto& trackTag : tracksTag) {
            if (clusterBaseParticleId != -1 && clusterBaseParticleId == trackTag.globalIndex()) {
              ptSum += clusterBase.energy() / std::cosh(clusterBase.eta());
              isClusterMatched = true;
              break;
            }
          }
          if (isClusterMatched) {
            break;
          }
        }
      }
    }
    if constexpr (jetsBaseIsMc) {
      for (const auto& trackBase : tracksBase) {
        if (std::find(particleTracker.begin(), particleTracker.end(), trackBase.globalIndex()) != particleTracker.end()) {
          continue;
        }
        auto trackBaseId = trackBase.globalIndex();
        for (const auto& clusterTag : clustersTag) {
          bool isClusterMatched = false;
          for (const auto& clusterTagParticleId : clusterTag.mcParticlesIds()) {
            if (trackBaseId != -1 && trackBaseId == clusterTagParticleId) {
              ptSum += trackBase.pt();
              isClusterMatched = true;
              break;
            }
          }
          if (isClusterMatched) {
            break;
          }
        }
      }
    }
  }
  if constexpr (isCandidate) {
    if constexpr (jetsTagIsMc) {
      for (auto const& candidateBase : candidatesBase) {
        if (jetcandidateutilities::isMatchedCandidate(candidateBase)) {
          const auto candidateBaseMcId = jetcandidateutilities::matchedParticleId(candidateBase, fullTracksBase, fullTracksTag);
          for (auto const& candidateTag : candidatesTag) {
            const auto candidateTagId = candidateTag.mcParticleId();
            if (candidateBaseMcId == candidateTagId) {
              ptSum += candidateBase.pt();
            }
            break; // should only be one
          }
        }
        break;
      }
    } else if constexpr (jetsBaseIsMc) {
      for (auto const& candidateTag : candidatesTag) {
        if (jetcandidateutilities::isMatchedCandidate(candidateTag)) {
          const auto candidateTagMcId = jetcandidateutilities::matchedParticleId(candidateTag, fullTracksTag, fullTracksBase);
          for (auto const& candidateBase : candidatesBase) {
            const auto candidateBaseId = candidateBase.mcParticleId();
            if (candidateTagMcId == candidateBaseId) {
              ptSum += candidateTag.pt();
            }
            break; // should only be one
          }
        }
        break;
      }
    } else {
      for (auto const& candidateBase : candidatesBase) {
        for (auto const& candidateTag : candidatesTag) {
          if (candidateBase.globalIndex() == candidateTag.globalIndex()) {
            ptSum += candidateBase.pt();
          }
          break; // should only be one
        }
        break;
      }
    }
  }
  return ptSum;
}

template <typename T, typename U>
auto getConstituents(T const& jet, U const& /*constituents*/)
{
  if constexpr (jetfindingutilities::isEMCALClusterTable<U>()) {
    return jet.template clusters_as<U>();
  } else if constexpr (jetcandidateutilities::isCandidateTable<U>() || jetcandidateutilities::isCandidateMcTable<U>()) {
    return jet.template candidates_as<U>();
  } else if constexpr (jetfindingutilities::isDummyTable<U>() || std::is_same_v<U, o2::aod::JCollisions> || std::is_same_v<U, o2::aod::JMcCollisions>) { // this is for the case where EMCal clusters or candidates are tested but no clusters or candidates exist and dummy tables are used, like in the case of charged jet analyses
    return nullptr;
  } else {
    return jet.template tracks_as<U>();
  }
}

template <bool jetsBaseIsMc, bool jetsTagIsMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q>
void MatchPt(T const& jetsBasePerCollision, U const& jetsTagPerCollision, std::vector<std::vector<int>>& baseToTagMatchingPt, std::vector<std::vector<int>>& tagToBaseMatchingPt, V const& tracksBase, M const& candidatesBase, N const& clustersBase, O const& tracksTag, P const& candidatesTag, Q const& clustersTag, float minPtFraction)
{
  float ptSumBase;
  float ptSumTag;
  for (const auto& jetBase : jetsBasePerCollision) {
    auto jetBaseTracks = getConstituents(jetBase, tracksBase);
    auto jetBaseClusters = getConstituents(jetBase, clustersBase);
    auto jetBaseCandidates = getConstituents(jetBase, candidatesBase);
    for (const auto& jetTag : jetsTagPerCollision) {
      if (std::round(jetBase.r()) != std::round(jetTag.r())) {
        continue;
      }
      auto jetTagTracks = getConstituents(jetTag, tracksTag);
      auto jetTagClusters = getConstituents(jetTag, clustersTag);
      auto jetTagCandidates = getConstituents(jetTag, candidatesTag);

      constexpr bool IsEMCAL{jetfindingutilities::isEMCALClusterTable<N>() || jetfindingutilities::isEMCALClusterTable<Q>()};
      constexpr bool IsCandidate{(jetcandidateutilities::isCandidateTable<M>() || jetcandidateutilities::isCandidateMcTable<M>()) && (jetcandidateutilities::isCandidateTable<P>() || jetcandidateutilities::isCandidateMcTable<P>())};
      ptSumBase = getPtSum<IsEMCAL, IsCandidate, jetsBaseIsMc, jetsTagIsMc>(jetBaseTracks, jetBaseCandidates, jetBaseClusters, jetTagTracks, jetTagCandidates, jetTagClusters, tracksBase, tracksTag);
      ptSumTag = getPtSum<IsEMCAL, IsCandidate, jetsTagIsMc, jetsBaseIsMc>(jetTagTracks, jetTagCandidates, jetTagClusters, jetBaseTracks, jetBaseCandidates, jetBaseClusters, tracksTag, tracksBase);
      if (ptSumBase > jetBase.pt() * minPtFraction) {
        baseToTagMatchingPt[jetBase.globalIndex()].push_back(jetTag.globalIndex());
      }
      if (ptSumTag > jetTag.pt() * minPtFraction) {
        tagToBaseMatchingPt[jetTag.globalIndex()].push_back(jetBase.globalIndex());
      }
    }
  }
}

// function that calls all the Match functions
template <bool jetsBaseIsMc, bool jetsTagIsMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename R>
void doAllMatching(T const& jetsBasePerCollision, U const& jetsTagPerCollision, std::vector<std::vector<int>>& baseToTagMatchingGeo, std::vector<std::vector<int>>& baseToTagMatchingPt, std::vector<std::vector<int>>& baseToTagMatchingHF, std::vector<std::vector<int>>& tagToBaseMatchingGeo, std::vector<std::vector<int>>& tagToBaseMatchingPt, std::vector<std::vector<int>>& tagToBaseMatchingHF, V const& candidatesBase, M const& tracksBase, N const& clustersBase, O const& candidatesTag, P const& tracksTag, R const& clustersTag, bool doMatchingGeo, bool doMatchingHf, bool doMatchingPt, float maxMatchingDistance, float minPtFraction)
{
  // geometric matching
  if (doMatchingGeo) {
    MatchGeo(jetsBasePerCollision, jetsTagPerCollision, baseToTagMatchingGeo, tagToBaseMatchingGeo, maxMatchingDistance);
  }
  // pt matching
  if (doMatchingPt) {
    MatchPt<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerCollision, jetsTagPerCollision, baseToTagMatchingPt, tagToBaseMatchingPt, tracksBase, candidatesBase, clustersBase, tracksTag, candidatesTag, clustersTag, minPtFraction);
  }
  // HF matching
  if constexpr (jetcandidateutilities::isCandidateTable<V>() || jetcandidateutilities::isCandidateMcTable<V>()) {
    if (doMatchingHf) {
      MatchHF<jetsBaseIsMc, jetsTagIsMc>(jetsBasePerCollision, jetsTagPerCollision, baseToTagMatchingHF, tagToBaseMatchingHF, candidatesBase, candidatesTag, tracksBase, tracksTag);
    }
  }
}

// function that does pair matching
template <bool jetsBaseIsMc, bool jetsTagIsMc, typename T, typename U, typename V, typename M, typename N, typename O>
void doPairMatching(T const& pairsBase, U const& pairsTag, std::vector<std::vector<int>>& baseToTagMatching, std::vector<std::vector<int>>& tagToBaseMatching, V const& /*candidatesBase*/, M const& tracksBase, N const& /*candidatesTag*/, O const& tracksTag)
{
  bool hasTrackBase1 = false;
  bool hasTrackBase2 = false;
  bool hasCandidateBase1 = false;
  bool hasCandidateBase2 = false;
  std::vector<int> pairsTagIndices;
  for (auto i = 0; i < pairsTag.size(); i++) {
    pairsTagIndices.push_back(i);
  }
  for (const auto& pairBase : pairsBase) {
    if (pairBase.has_track1()) {
      hasTrackBase1 = true;
    }
    if (pairBase.has_track2()) {
      hasTrackBase2 = true;
    }
    if (pairBase.has_candidate1()) {
      hasCandidateBase1 = true;
    }
    if (pairBase.has_candidate2()) {
      hasCandidateBase2 = true;
    }
    int matchedPairTagIndex = -1;
    for (auto pairTagIndex : pairsTagIndices) {
      const auto& pairTag = pairsTag.iteratorAt(pairTagIndex);
      if (hasTrackBase1 && !pairTag.has_track1()) {
        continue;
      }
      if (hasTrackBase2 && !pairTag.has_track2()) {
        continue;
      }
      if (hasCandidateBase1 && !pairTag.has_candidate1()) {
        continue;
      }
      if (hasCandidateBase2 && !pairTag.has_candidate2()) {
        continue;
      }
      int nMatched = 0;
      bool isMatched = false;
      if (hasTrackBase1) {
        const auto& trackBase1 = pairBase.template track1_as<M>();
        const auto& trackTag1 = pairTag.template track1_as<O>();
        if constexpr (jetsTagIsMc) {
          if (trackBase1.mcParticleId() == trackTag1.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        } else if constexpr (jetsBaseIsMc) {
          if (trackBase1.globalIndex() == trackTag1.mcParticleId()) {
            nMatched++;
            isMatched = true;
          }
        } else {
          if (trackBase1.globalIndex() == trackTag1.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        }
        if (!isMatched) {
          continue;
        }
      }
      isMatched = false;

      if (hasTrackBase2) {
        const auto& trackBase2 = pairBase.template track2_as<M>();
        const auto& trackTag2 = pairTag.template track2_as<O>();
        if constexpr (jetsTagIsMc) {
          if (trackBase2.mcParticleId() == trackTag2.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        } else if constexpr (jetsBaseIsMc) {
          if (trackBase2.globalIndex() == trackTag2.mcParticleId()) {
            nMatched++;
            isMatched = true;
          }
        } else {
          if (trackBase2.globalIndex() == trackTag2.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        }
        if (!isMatched) {
          continue;
        }
      }
      isMatched = false;
      if (hasCandidateBase1) {
        const auto& candidateBase1 = pairBase.template candidate1_as<V>();
        const auto& candidateTag1 = pairTag.template candidate1_as<N>();
        if constexpr (jetsTagIsMc) {
          if (jetcandidateutilities::isMatchedCandidate(candidateBase1)) {
            const auto candidateBaseMcId = jetcandidateutilities::matchedParticleId(candidateBase1, tracksBase, tracksTag);
            if (candidateBaseMcId == candidateTag1.globalIndex()) {
              nMatched++;
              isMatched = true;
            }
          }
        } else if constexpr (jetsBaseIsMc) {
          if (jetcandidateutilities::isMatchedCandidate(candidateTag1)) {
            const auto candidateTagMcId = jetcandidateutilities::matchedParticleId(candidateTag1, tracksTag, tracksBase);
            if (candidateTagMcId == candidateBase1.globalIndex()) {
              nMatched++;
              isMatched = true;
            }
          }
        } else {
          if (candidateBase1.globalIndex() == candidateTag1.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        }
        if (!isMatched) {
          continue;
        }
      }
      isMatched = false;
      if (hasCandidateBase2) {
        const auto& candidateBase2 = pairBase.template candidate2_as<V>();
        const auto& candidateTag2 = pairTag.template candidate2_as<N>();
        if constexpr (jetsTagIsMc) {
          if (jetcandidateutilities::isMatchedCandidate(candidateBase2)) {
            const auto candidateBaseMcId = jetcandidateutilities::matchedParticleId(candidateBase2, tracksBase, tracksTag);
            if (candidateBaseMcId == candidateTag2.globalIndex()) {
              nMatched++;
              isMatched = true;
            }
          }
        } else if constexpr (jetsBaseIsMc) {
          if (jetcandidateutilities::isMatchedCandidate(candidateTag2)) {
            const auto candidateTagMcId = jetcandidateutilities::matchedParticleId(candidateTag2, tracksTag, tracksBase);
            if (candidateTagMcId == candidateBase2.globalIndex()) {
              nMatched++;
              isMatched = true;
            }
          }
        } else {
          if (candidateBase2.globalIndex() == candidateTag2.globalIndex()) {
            nMatched++;
            isMatched = true;
          }
        }
        if (!isMatched) {
          continue;
        }
      }

      if (nMatched == 2) {
        baseToTagMatching[pairBase.globalIndex()].push_back(pairTag.globalIndex());
        tagToBaseMatching[pairTag.globalIndex()].push_back(pairBase.globalIndex());
        matchedPairTagIndex = pairTagIndex;
        break; // can only be one match per jet
      }
    }
    if (matchedPairTagIndex != -1) {
      pairsTagIndices.erase(std::find(pairsTagIndices.begin(), pairsTagIndices.end(), matchedPairTagIndex));
    }
  }
}

}; // namespace jetmatchingutilities
#endif // PWGJE_CORE_JETMATCHINGUTILITIES_H_
