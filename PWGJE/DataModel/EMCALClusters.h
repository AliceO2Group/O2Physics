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
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL

#ifndef PWGJE_DATAMODEL_EMCALCLUSTERS_H_
#define PWGJE_DATAMODEL_EMCALCLUSTERS_H_

#include <string>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "EMCALClusterDefinition.h"

namespace o2::aod
{
namespace emcalcluster
{

// define global cluster definitions
// New definitions should be added here!
const EMCALClusterDefinition kV3NoSplit(ClusterAlgorithm_t::kV3, 0, 1, "kV3NoSplit", 0.5, 0.1, -10000, 10000, 20000, false, 0., false);
const EMCALClusterDefinition kV3NoSplitLowSeed(ClusterAlgorithm_t::kV3, 1, 1, "kV3NoSplitLowSeed", 0.3, 0.1, -10000, 10000, 20000, false, 0., false);
const EMCALClusterDefinition kV3NoSplitLowerSeed(ClusterAlgorithm_t::kV3, 2, 1, "kV3NoSplitLowerSeed", 0.2, 0.1, -10000, 10000, 20000, false, 0., false);
const EMCALClusterDefinition kV3Default(ClusterAlgorithm_t::kV3, 10, 1, "kV3Default", 0.5, 0.1, -10000, 10000, 20000, true, 0.03, false);
const EMCALClusterDefinition kV3MostSplit(ClusterAlgorithm_t::kV3, 11, 1, "kV3MostSplit", 0.5, 0.1, -10000, 10000, 20000, true, 0., false);
const EMCALClusterDefinition kV3LowSeed(ClusterAlgorithm_t::kV3, 12, 1, "kV3LowSeed", 0.3, 0.1, -10000, 10000, 20000, true, 0.03, false);
const EMCALClusterDefinition kV3MostSplitLowSeed(ClusterAlgorithm_t::kV3, 13, 1, "kV3MostSplitLowSeed", 0.3, 0.1, -10000, 10000, 20000, true, 0., false);
const EMCALClusterDefinition kV3StrictTime(ClusterAlgorithm_t::kV3, 20, 1, "kV3StrictTime", 0.5, 0.1, -500, 500, 20000, true, 0.03, false);
const EMCALClusterDefinition kV3StricterTime(ClusterAlgorithm_t::kV3, 21, 1, "kV3StricterTime", 0.5, 0.1, -100, 100, 20000, true, 0.03, false);
const EMCALClusterDefinition kV3MostStrictTime(ClusterAlgorithm_t::kV3, 22, 1, "kV3MostStrictTime", 0.5, 0.1, -50, 50, 20000, true, 0.03, false);
const EMCALClusterDefinition kV3Default5x5(ClusterAlgorithm_t::kV3, 30, 1, "kV3Default5x5", 0.5, 0.1, -10000, 10000, 20000, true, 0.03, true);
const EMCALClusterDefinition kV3SmallTimeDiff(ClusterAlgorithm_t::kV3, 40, 1, "kV3SmallTimeDiff", 0.5, 0.1, -10000, 10000, 500, true, 0.03, false);
const EMCALClusterDefinition kV3SmallerTimeDiff(ClusterAlgorithm_t::kV3, 41, 1, "kV3SmallerTimeDiff", 0.5, 0.1, -10000, 10000, 100, true, 0.03, false);
const EMCALClusterDefinition kV3SmallestTimeDiff(ClusterAlgorithm_t::kV3, 42, 1, "kV3SmallestTimeDiff", 0.5, 0.1, -10000, 10000, 50, true, 0.03, false);
const EMCALClusterDefinition kV3MostSplitSmallTimeDiff(ClusterAlgorithm_t::kV3, 43, 1, "kV3MostSplitSmallTimeDiff", 0.5, 0.1, -10000, 10000, 500, true, 0., false);
const EMCALClusterDefinition kV3MostSplitSmallerTimeDiff(ClusterAlgorithm_t::kV3, 44, 1, "kV3MostSplitSmallerTimeDiff", 0.5, 0.1, -10000, 10000, 100, true, 0., false);
const EMCALClusterDefinition kV3MostSplitSmallestTimeDiff(ClusterAlgorithm_t::kV3, 45, 1, "kV3MostSplitSmallestTimeDiff", 0.5, 0.1, -10000, 10000, 50, true, 0., false);

/// \brief function returns EMCALClusterDefinition for the given name
/// \param name name of the cluster definition
/// \return EMCALClusterDefinition for the given name
const EMCALClusterDefinition getClusterDefinitionFromString(const std::string& clusterDefinitionName)
{
  if (clusterDefinitionName == "kV3NoSplit") {
    return kV3NoSplit;
  } else if (clusterDefinitionName == "kV3NoSplitLowSeed") {
    return kV3NoSplitLowSeed;
  } else if (clusterDefinitionName == "kV3NoSplitLowerSeed") {
    return kV3NoSplitLowerSeed;
  } else if (clusterDefinitionName == "kV3Default") {
    return kV3Default;
  } else if (clusterDefinitionName == "kV3MostSplit") {
    return kV3MostSplit;
  } else if (clusterDefinitionName == "kV3LowSeed") {
    return kV3LowSeed;
  } else if (clusterDefinitionName == "kV3MostSplitLowSeed") {
    return kV3MostSplitLowSeed;
  } else if (clusterDefinitionName == "kV3StrictTime") {
    return kV3StrictTime;
  } else if (clusterDefinitionName == "kV3StricterTime") {
    return kV3StricterTime;
  } else if (clusterDefinitionName == "kV3MostStrictTime") {
    return kV3MostStrictTime;
  } else if (clusterDefinitionName == "kV3Default5x5") {
    return kV3Default5x5;
  } else if (clusterDefinitionName == "kV3SmallTimeDiff") {
    return kV3SmallTimeDiff;
  } else if (clusterDefinitionName == "kV3SmallerTimeDiff") {
    return kV3SmallerTimeDiff;
  } else if (clusterDefinitionName == "kV3SmallestTimeDiff") {
    return kV3SmallestTimeDiff;
  } else if (clusterDefinitionName == "kV3MostSplitSmallTimeDiff") {
    return kV3MostSplitSmallTimeDiff;
  } else if (clusterDefinitionName == "kV3MostSplitSmallerTimeDiff") {
    return kV3MostSplitSmallerTimeDiff;
  } else if (clusterDefinitionName == "kV3MostSplitSmallestTimeDiff") {
    return kV3MostSplitSmallestTimeDiff;
  } else {
    throw std::invalid_argument("Cluster definition name not recognized");
  }
};

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
DECLARE_SOA_INDEX_COLUMN(Track, track); //! linked to Track table only for tracks that were matched
} // namespace emcalmatchedtrack
DECLARE_SOA_TABLE(EMCALMatchedTracks, "AOD", "EMCMATCHTRACKS",                                     //!
                  o2::soa::Index<>, emcalclustercell::EMCALClusterId, emcalmatchedtrack::TrackId); //!
using EMCALMatchedTrack = EMCALMatchedTracks::iterator;
} // namespace o2::aod
#endif // PWGJE_DATAMODEL_EMCALCLUSTERS_H_
