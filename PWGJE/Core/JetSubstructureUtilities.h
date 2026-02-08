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

/// \file JetSubstructureUtilities.h
/// \brief Jet Substructure related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETSUBSTRUCTUREUTILITIES_H_
#define PWGJE_CORE_JETSUBSTRUCTUREUTILITIES_H_

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetCandidateUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

#include <Framework/ASoA.h>

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/MeasureDefinition.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <cmath>
#include <type_traits>
#include <vector>

namespace jetsubstructureutilities
{

/**
 * convert an O2Physics jet to a fastjet pseudojet object, returning its clusterSequence
 *
 * @param jet jet to be converted
 * @param tracks vector of constituent tracks
 * @param clusters vector of constituent clusters
 * @param candidates vector of constituent candidates
 * @param pseudoJet converted pseudoJet object which is passed by reference
 */
template <typename T, typename U, typename V, typename O>
fastjet::ClusterSequenceArea jetToPseudoJet(T const& jet, U const& /*tracks*/, V const& /*clusters*/, O const& /*candidates*/, fastjet::PseudoJet& pseudoJet, int hadronicCorrectionType = 0)
{
  std::vector<fastjet::PseudoJet> jetConstituents;
  for (auto& jetConstituent : jet.template tracks_as<U>()) {
    fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
  }
  if constexpr (std::is_same_v<std::decay_t<typename V::iterator>, o2::aod::JetClusters::iterator> || std::is_same_v<std::decay_t<typename V::iterator>, o2::aod::JetClusters::filtered_iterator>) {
    for (auto& jetClusterConstituent : jet.template clusters_as<V>()) {
      fastjetutilities::fillClusters(jetClusterConstituent, jetConstituents, jetClusterConstituent.globalIndex(), hadronicCorrectionType);
    }
  }
  if constexpr (jetcandidateutilities::isCandidateTable<O>() || jetcandidateutilities::isCandidateMcTable<O>()) {
    for (auto& jetHFConstituent : jet.template candidates_as<O>()) {
      fastjetutilities::fillTracks(jetHFConstituent, jetConstituents, jetHFConstituent.globalIndex(), JetConstituentStatus::candidate, jetcandidateutilities::getTablePDGMass<O>());
    }
  }
  std::vector<fastjet::PseudoJet> jetReclustered;

  JetFinder jetReclusterer;
  jetReclusterer.isReclustering = true;
  jetReclusterer.jetR = jet.r() / 100.0;
  fastjet::ClusterSequenceArea clusterSeq = jetReclusterer.findJets(jetConstituents, jetReclustered);
  jetReclustered = sorted_by_pt(jetReclustered);
  pseudoJet = jetReclustered[0];
  return clusterSeq;
}

/**
 * returns a vector with Nsubjettiness variables
 *
 * @param jet jet
 * @param tracks track table to be added
 * @param clusters clusters table to be added (if no clusters just add track table here)
 * @param candidates candidates table to be added (if no candidates just add track table here)
 * @param nMax returns a vector filled with TauN values upto N (the first entry is the distance between axes in tau2)
 * @param reclusteringAlgorithm type of reclustering algorithm used to find Nsubjettiness axes
 * @param doSoftDrop apply SoftDrop
 * @param zCut minimim momentum sharing fraction needed to satisfy the SoftDrop condition
 * @param beta angular exponent in the SoftDrop condition
 */

// function that returns the N-subjettiness ratio and the distance betewwen the two axes considered for tau2, in the form of a vector
template <typename T, typename U, typename V, typename O, typename M>
std::vector<float> getNSubjettiness(T const& jet, U const& tracks, V const& clusters, O const& candidates, std::vector<fastjet::PseudoJet>::size_type nMax, M const& reclusteringAlgorithm, bool doSoftDrop = false, float zCut = 0.1, float beta = 0.0, int hadronicCorrectionType = 0)
{
  std::vector<float> result;
  fastjet::PseudoJet pseudoJet;
  fastjet::ClusterSequenceArea clusterSeq(jetToPseudoJet(jet, tracks, clusters, candidates, pseudoJet, hadronicCorrectionType));
  if (doSoftDrop) {
    fastjet::contrib::SoftDrop softDrop(beta, zCut);
    pseudoJet = softDrop(pseudoJet);
  }

  for (std::vector<fastjet::PseudoJet>::size_type n = 0; n < nMax + 1; n++) {
    result.push_back(-1.0 * (n + 1));
  }

  for (std::vector<fastjet::PseudoJet>::size_type n = 1; n <= nMax; n++) {
    if (pseudoJet.constituents().size() < n) { // Tau_N needs at least N tracks
      return result;
    }
    fastjet::contrib::Nsubjettiness nSub(n, reclusteringAlgorithm, fastjet::contrib::NormalizedMeasure(1.0, jet.r() / 100.0));
    result[n] = nSub.result(pseudoJet);
    if (n == 2) {
      std::vector<fastjet::PseudoJet> nSubAxes = nSub.currentAxes(); // gets the two axes used in the 2-subjettiness calculation
      result[0] = nSubAxes[0].delta_R(nSubAxes[1]);                  // distance between axes for 2-subjettiness
    }
  }
  return result;
}

}; // namespace jetsubstructureutilities

#endif // PWGJE_CORE_JETSUBSTRUCTUREUTILITIES_H_
