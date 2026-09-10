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

/// \file EMCALClusters.h
/// \brief Table definitions for EMCAL analysis clusters
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL, Florian Jonas <florian.jonas@cern.ch>, Marvin Hemmer <marvin.hemmer@cern.ch>

#ifndef PWGJE_DATAMODEL_EMCALCLUSTERS_H_
#define PWGJE_DATAMODEL_EMCALCLUSTERS_H_

#include "EMCALClusterDefinition.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::aod
{
namespace emcalcluster
{

// define global cluster definitions
// New definitions should be added here!
inline const EMCALClusterDefinition kV3NoSplit(ClusterAlgorithm::kV3, 0, 1, "kV3NoSplit", 0.5, 0.1, -10000, 10000, 20000, false, 0., false);
inline const EMCALClusterDefinition kV3NoSplitLowSeed(ClusterAlgorithm::kV3, 1, 1, "kV3NoSplitLowSeed", 0.3, 0.1, -10000, 10000, 20000, false, 0., false);
inline const EMCALClusterDefinition kV3NoSplitLowerSeed(ClusterAlgorithm::kV3, 2, 1, "kV3NoSplitLowerSeed", 0.2, 0.1, -10000, 10000, 20000, false, 0., false);
inline const EMCALClusterDefinition kV3Default(ClusterAlgorithm::kV3, 10, 1, "kV3Default", 0.5, 0.1, -10000, 10000, 20000, true, 0.03, false);
inline const EMCALClusterDefinition kV3MostSplit(ClusterAlgorithm::kV3, 11, 1, "kV3MostSplit", 0.5, 0.1, -10000, 10000, 20000, true, 0., false);
inline const EMCALClusterDefinition kV3LowSeed(ClusterAlgorithm::kV3, 12, 1, "kV3LowSeed", 0.3, 0.1, -10000, 10000, 20000, true, 0.03, false);
inline const EMCALClusterDefinition kV3MostSplitLowSeed(ClusterAlgorithm::kV3, 13, 1, "kV3MostSplitLowSeed", 0.3, 0.1, -10000, 10000, 20000, true, 0., false);
inline const EMCALClusterDefinition kV3StrictTime(ClusterAlgorithm::kV3, 20, 1, "kV3StrictTime", 0.5, 0.1, -500, 500, 20000, true, 0.03, false);
inline const EMCALClusterDefinition kV3StricterTime(ClusterAlgorithm::kV3, 21, 1, "kV3StricterTime", 0.5, 0.1, -100, 100, 20000, true, 0.03, false);
inline const EMCALClusterDefinition kV3MostStrictTime(ClusterAlgorithm::kV3, 22, 1, "kV3MostStrictTime", 0.5, 0.1, -50, 50, 20000, true, 0.03, false);
inline const EMCALClusterDefinition kV3Default5x5(ClusterAlgorithm::kV3, 30, 1, "kV3Default5x5", 0.5, 0.1, -10000, 10000, 20000, true, 0.03, true);
inline const EMCALClusterDefinition kV3SmallTimeDiff(ClusterAlgorithm::kV3, 40, 1, "kV3SmallTimeDiff", 0.5, 0.1, -10000, 10000, 500, true, 0.03, false);
inline const EMCALClusterDefinition kV3SmallerTimeDiff(ClusterAlgorithm::kV3, 41, 1, "kV3SmallerTimeDiff", 0.5, 0.1, -10000, 10000, 100, true, 0.03, false);
inline const EMCALClusterDefinition kV3SmallestTimeDiff(ClusterAlgorithm::kV3, 42, 1, "kV3SmallestTimeDiff", 0.5, 0.1, -10000, 10000, 50, true, 0.03, false);
inline const EMCALClusterDefinition kV3MostSplitSmallTimeDiff(ClusterAlgorithm::kV3, 43, 1, "kV3MostSplitSmallTimeDiff", 0.5, 0.1, -10000, 10000, 500, true, 0., false);
inline const EMCALClusterDefinition kV3MostSplitSmallerTimeDiff(ClusterAlgorithm::kV3, 44, 1, "kV3MostSplitSmallerTimeDiff", 0.5, 0.1, -10000, 10000, 100, true, 0., false);
inline const EMCALClusterDefinition kV3MostSplitSmallestTimeDiff(ClusterAlgorithm::kV3, 45, 1, "kV3MostSplitSmallestTimeDiff", 0.5, 0.1, -10000, 10000, 50, true, 0., false);
inline const EMCALClusterDefinition kV3MostSplitSmallestTimeDiffLowestSeed(ClusterAlgorithm::kV3, 50, 1, "kV3MostSplitSmallestTimeDiffLowestSeed", 0.1, 0.1, -10000, 10000, 50, true, 0., false);
inline const EMCALClusterDefinition kV3MostSplitSmallestTimeDiffLowSeed(ClusterAlgorithm::kV3, 51, 1, "kV3MostSplitSmallestTimeDiffLowSeed", 0.3, 0.1, -10000, 10000, 50, true, 0., false);
inline const EMCALClusterDefinition kV3MostSplitSmallestTimeDiffLowerSeed(ClusterAlgorithm::kV3, 52, 1, "kV3MostSplitSmallestTimeDiffLowerSeed", 0.2, 0.1, -10000, 10000, 50, true, 0., false);

/// \brief function returns EMCALClusterDefinition for the given storage ID
/// \param storageID storage ID of the cluster definition
/// \return EMCALClusterDefinition for the given storage ID
inline const EMCALClusterDefinition& getClusterDefinitionFromID(int storageID)
{
  switch (storageID) {
    case 0:
      return kV3NoSplit;
    case 1:
      return kV3NoSplitLowSeed;
    case 2:
      return kV3NoSplitLowerSeed;
    case 10:
      return kV3Default;
    case 11:
      return kV3MostSplit;
    case 12:
      return kV3LowSeed;
    case 13:
      return kV3MostSplitLowSeed;
    case 20:
      return kV3StrictTime;
    case 21:
      return kV3StricterTime;
    case 22:
      return kV3MostStrictTime;
    case 30:
      return kV3Default5x5;
    case 40:
      return kV3SmallTimeDiff;
    case 41:
      return kV3SmallerTimeDiff;
    case 42:
      return kV3SmallestTimeDiff;
    case 43:
      return kV3MostSplitSmallTimeDiff;
    case 44:
      return kV3MostSplitSmallerTimeDiff;
    case 45:
      return kV3MostSplitSmallestTimeDiff;
    case 50:
      return kV3MostSplitSmallestTimeDiffLowestSeed;
    case 51:
      return kV3MostSplitSmallestTimeDiffLowSeed;
    case 52:
      return kV3MostSplitSmallestTimeDiffLowerSeed;
    default:
      throw std::invalid_argument("Cluster definition storageID not recognized: " + std::to_string(storageID));
  }
}

/// \brief function returns EMCALClusterDefinition for the given name
/// \param clusterDefinitionName name of the cluster definition
/// \return EMCALClusterDefinition for the given name
inline const EMCALClusterDefinition& getClusterDefinitionFromString(const std::string& clusterDefinitionName)
{
  static const std::unordered_map<std::string, int> nameToID = {
    {"kV3NoSplit", 0},
    {"kV3NoSplitLowSeed", 1},
    {"kV3NoSplitLowerSeed", 2},
    {"kV3Default", 10},
    {"kV3MostSplit", 11},
    {"kV3LowSeed", 12},
    {"kV3MostSplitLowSeed", 13},
    {"kV3StrictTime", 20},
    {"kV3StricterTime", 21},
    {"kV3MostStrictTime", 22},
    {"kV3Default5x5", 30},
    {"kV3SmallTimeDiff", 40},
    {"kV3SmallerTimeDiff", 41},
    {"kV3SmallestTimeDiff", 42},
    {"kV3MostSplitSmallTimeDiff", 43},
    {"kV3MostSplitSmallerTimeDiff", 44},
    {"kV3MostSplitSmallestTimeDiff", 45},
    {"kV3MostSplitSmallestTimeDiffLowestSeed", 50},
    {"kV3MostSplitSmallestTimeDiffLowSeed", 51},
    {"kV3MostSplitSmallestTimeDiffLowerSeed", 52},
  };

  auto it = nameToID.find(clusterDefinitionName);
  if (it == nameToID.end()) {
    throw std::invalid_argument("Cluster definition name not recognized: " + clusterDefinitionName);
  }
  return getClusterDefinitionFromID(it->second);
}

DECLARE_SOA_INDEX_COLUMN(Collision, collision);                        //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);                                      //! bunch crossing ID used as index for ambiguous clusters
DECLARE_SOA_COLUMN(ID, id, int);                                       //! cluster ID identifying cluster in event
DECLARE_SOA_COLUMN(Energy, energy, float);                             //! cluster energy (GeV)
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                     //! cluster core energy (GeV)
DECLARE_SOA_COLUMN(RawEnergy, rawEnergy, float);                       //! raw cluster energy (GeV)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                                   //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                   //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, int);                               //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                 //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                          //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //! distance to bad channel
DECLARE_SOA_COLUMN(NLM, nlm, int);                                     //! number of local maxima
DECLARE_SOA_COLUMN(Definition, definition, int);                       //! cluster definition, see EMCALClusterDefinition.h

} // namespace emcalcluster
// table of clusters that could be matched to a collision
DECLARE_SOA_TABLE(EMCALClusters, "AOD", "EMCALCLUSTERS", //!
                  o2::soa::Index<>, emcalcluster::CollisionId, emcalcluster::ID, emcalcluster::Energy,
                  emcalcluster::CoreEnergy, emcalcluster::RawEnergy, emcalcluster::Eta, emcalcluster::Phi,
                  emcalcluster::M02, emcalcluster::M20, emcalcluster::NCells, emcalcluster::Time,
                  emcalcluster::IsExotic, emcalcluster::DistanceToBadChannel, emcalcluster::NLM, emcalcluster::Definition);
// table of ambiguous clusters that could not be matched to a collision
DECLARE_SOA_TABLE(EMCALAmbiguousClusters, "AOD", "EMCALAMBCLUS", //!
                  o2::soa::Index<>, emcalcluster::BCId, emcalcluster::ID, emcalcluster::Energy,
                  emcalcluster::CoreEnergy, emcalcluster::RawEnergy, emcalcluster::Eta, emcalcluster::Phi,
                  emcalcluster::M02, emcalcluster::M20, emcalcluster::NCells, emcalcluster::Time,
                  emcalcluster::IsExotic, emcalcluster::DistanceToBadChannel, emcalcluster::NLM, emcalcluster::Definition);

using EMCALCluster = EMCALClusters::iterator;
using EMCALAmbiguousCluster = EMCALAmbiguousClusters::iterator;

namespace emcalclustermc
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(McParticle, mcParticle);         //! Array of MC particles that deposited energy in this calo cell
DECLARE_SOA_COLUMN(AmplitudeA, amplitudeA, std::vector<float>); //! Energy fraction deposited by a particle inside this calo cell.
} // namespace emcalclustermc
// table of cluster MC info that could be matched to a collision
DECLARE_SOA_TABLE(EMCALMCClusters, "AOD", "EMCALMCCLUSTERS", //!
                  emcalclustermc::McParticleIds, emcalclustermc::AmplitudeA);

using EMCALMCCluster = EMCALMCClusters::iterator;

// table of cluster MC info that could not be matched to a collision
DECLARE_SOA_TABLE(EMCALAmbiguousMCClusters, "AOD", "EMCALAMBMCCLS", //!
                  emcalclustermc::McParticleIds, emcalclustermc::AmplitudeA);

using EMCALAmbiguousMCCluster = EMCALAmbiguousMCClusters::iterator;

namespace emcalclustercell
{
// declare index column pointing to cluster table
DECLARE_SOA_INDEX_COLUMN(EMCALCluster, emcalcluster); //! linked to EMCalClusters table
DECLARE_SOA_INDEX_COLUMN(Calo, calo);                 //! linked to calo cells

// declare index column pointing to ambiguous cluster table
DECLARE_SOA_INDEX_COLUMN(EMCALAmbiguousCluster, emcalambiguouscluster); //! linked to EMCalAmbiguousClusters table
} // namespace emcalclustercell
DECLARE_SOA_TABLE(EMCALClusterCells, "AOD", "EMCCLUSCELLS",                                               //!
                  o2::soa::Index<>, emcalclustercell::EMCALClusterId, emcalclustercell::CaloId);          //!
DECLARE_SOA_TABLE(EMCALAmbiguousClusterCells, "AOD", "EMCAMBBCLUSCLS",                                    //!
                  o2::soa::Index<>, emcalclustercell::EMCALAmbiguousClusterId, emcalclustercell::CaloId); //!
using EMCALClusterCell = EMCALClusterCells::iterator;
using EMCALAmbiguousClusterCell = EMCALAmbiguousClusterCells::iterator;
namespace emcalmatchedtrack
{
DECLARE_SOA_INDEX_COLUMN(Track, track);        //! linked to Track table only for tracks that were matched
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float); //! difference between matched track and cluster azimuthal angle
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float); //! difference between matched track and cluster pseudorapidity
} // namespace emcalmatchedtrack
DECLARE_SOA_TABLE(EMCALMatchedTracks, "AOD", "EMCMATCHTRACKS", //!
                  o2::soa::Index<>, emcalclustercell::EMCALClusterId, emcalmatchedtrack::TrackId,
                  emcalmatchedtrack::DeltaPhi, emcalmatchedtrack::DeltaEta); //!
using EMCALMatchedTrack = EMCALMatchedTracks::iterator;

// table for matched secondary tracks
DECLARE_SOA_TABLE(EMCMatchSecs, "AOD", "EMCMATCHSEC", //!
                  o2::soa::Index<>, emcalclustercell::EMCALClusterId, emcalmatchedtrack::TrackId,
                  emcalmatchedtrack::DeltaPhi, emcalmatchedtrack::DeltaEta); //!
using EMCMatchSec = EMCMatchSecs::iterator;
} // namespace o2::aod
#endif // PWGJE_DATAMODEL_EMCALCLUSTERS_H_
