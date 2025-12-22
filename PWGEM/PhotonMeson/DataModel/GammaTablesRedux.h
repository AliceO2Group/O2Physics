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
///
/// \file GammaTablesRedux.h
/// \brief This header provides the table definitions to store very lightweight EMCal clusters per collision
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) - Goethe University Frankfurt
///

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLESREDUX_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLESREDUX_H_

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <fairlogger/Logger.h>

#include <cmath>
#include <cstdint>
#include <limits>

namespace o2::aod
{

namespace emcdownscaling
{
enum Observable {
  kDefinition,
  kEnergy,
  kEta,
  kPhi,
  kNCellsExo,
  kM02,
  kTime,
  kDeltaEta,
  kDeltaPhi,
  nObservables
};

// Values in tables are stored in downscaled format to save disk space
const float downscalingFactors[nObservables]{
  1E0, // Cluster definition
  1E3, // Cluster energy
  1E4, // Cluster eta
  1E4, // Cluster phi
  1E0, // Number of cells + exotic
  1E4, // M02
  1E2, // Cluster time
  1E5, // diff between cluster and track in eta
  1E5, // diff between cluster and track in phi
};

/// \brief convert values for storage into correspoding smaller types
/// \param valueIn input value that should be stored
/// \param observable type of observable that should be stored to determine downscaling
template <typename OutputType, typename InputType>
OutputType convertForStorage(InputType const& valueIn, Observable observable)
{
  InputType valueToBeChecked = std::lround(valueIn * downscalingFactors[observable]);
  if (valueToBeChecked < std::numeric_limits<OutputType>::lowest()) {
    LOG(warning) << "Value " << valueToBeChecked << " of observable " << observable << " below lowest possible value of " << typeid(OutputType).name() << ": " << static_cast<float>(std::numeric_limits<OutputType>::lowest());
    valueToBeChecked = static_cast<float>(std::numeric_limits<OutputType>::lowest());
  }
  if (valueToBeChecked > std::numeric_limits<OutputType>::max()) {
    LOG(warning) << "Value " << valueToBeChecked << " of observable " << observable << " obove highest possible value of " << typeid(OutputType).name() << ": " << static_cast<float>(std::numeric_limits<OutputType>::max());
    valueToBeChecked = static_cast<float>(std::numeric_limits<OutputType>::max());
  }
  return static_cast<OutputType>(valueToBeChecked);
}
} // namespace emcdownscaling

namespace mincluster
{
DECLARE_SOA_COLUMN(StoredDefinition, storedDefinition, int8_t); //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_COLUMN(StoredE, storedE, uint16_t);                 //! cluster energy (1 MeV -> Maximum cluster energy of ~65 GeV)
DECLARE_SOA_COLUMN(StoredP, storedP, uint16_t);                 //! track momentum (1 MeV -> Maximum track momentum of ~65 GeV)
DECLARE_SOA_COLUMN(StoredPt, storedPt, uint16_t);               //! track pT (1 MeV -> Maximum track pT of ~65 GeV)
DECLARE_SOA_COLUMN(StoredEta, storedEta, int16_t);              //! cluster pseudorapidity (x10,000)
DECLARE_SOA_COLUMN(StoredPhi, storedPhi, uint16_t);             //! cluster azimuthal angle (x10 000) from 0 to 2pi
DECLARE_SOA_COLUMN(StoredNCells, storedNCells, uint8_t);        //! number of cells in cluster last 7 digits, first digit for exotic flag
DECLARE_SOA_COLUMN(StoredM02, storedM02, int16_t);              //! shower shape long axis (x10 000)
DECLARE_SOA_COLUMN(StoredTime, storedTime, int16_t);            //! cluster time (10 ps resolution)

DECLARE_SOA_COLUMN(StoredDeltaEta, storedDeltaEta, int16_t); //! cluster track difference pseudorapidity (x100,000)
DECLARE_SOA_COLUMN(StoredDeltaPhi, storedDeltaPhi, int16_t); //! cluster track difference azimuthal angle (x100,000)

DECLARE_SOA_DYNAMIC_COLUMN(Definition, definition, [](int8_t storedDefinition) -> int8_t { return storedDefinition; });                                          //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](uint16_t storedE) -> float { return storedE / emcdownscaling::downscalingFactors[emcdownscaling::kEnergy]; });               //! cluster energy (GeV)
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](int16_t storedEta) -> float { return storedEta / emcdownscaling::downscalingFactors[emcdownscaling::kEta]; });           //! cluster pseudorapidity
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](uint16_t storedPhi) -> float { return storedPhi / emcdownscaling::downscalingFactors[emcdownscaling::kPhi]; });          //! cluster azimuthal angle (0 to 2pi)
DECLARE_SOA_DYNAMIC_COLUMN(NCells, nCells, [](uint8_t storedNCells) -> uint8_t { return storedNCells & 0x7F; });                                                 //! number of cells in cluster
DECLARE_SOA_DYNAMIC_COLUMN(M02, m02, [](int16_t storedM02) -> float { return storedM02 / emcdownscaling::downscalingFactors[emcdownscaling::kM02]; });           //! shower shape long axis
DECLARE_SOA_DYNAMIC_COLUMN(Time, time, [](int16_t storedTime) -> float { return storedTime / emcdownscaling::downscalingFactors[emcdownscaling::kTime]; });      //! cluster time (ns)
DECLARE_SOA_DYNAMIC_COLUMN(IsExotic, isExotic, [](int16_t storedNCells) -> bool { return static_cast<bool>(storedNCells & 0x80); });                             //! flag to mark cluster as exotic
DECLARE_SOA_DYNAMIC_COLUMN(TrackP, trackP, [](uint16_t storedP) -> float { return storedP / emcdownscaling::downscalingFactors[emcdownscaling::kEnergy]; });     //! track momentum (GeV)
DECLARE_SOA_DYNAMIC_COLUMN(TrackPt, trackPt, [](uint16_t storedPt) -> float { return storedPt / emcdownscaling::downscalingFactors[emcdownscaling::kEnergy]; }); //! track pT (GeV)

DECLARE_SOA_DYNAMIC_COLUMN(DeltaEta, deltaEta, [](int16_t storedDeltaEta) -> float { return storedDeltaEta / emcdownscaling::downscalingFactors[emcdownscaling::kDeltaEta]; }); //! cluster track difference pseudorapidity
DECLARE_SOA_DYNAMIC_COLUMN(DeltaPhi, deltaPhi, [](int16_t storedDeltaPhi) -> float { return storedDeltaPhi / emcdownscaling::downscalingFactors[emcdownscaling::kDeltaPhi]; }); //! cluster track difference azimuthal angle

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float storedE, float storedEta) -> float { return storedE / emcdownscaling::downscalingFactors[emcdownscaling::kEnergy] / std::cosh(storedEta / emcdownscaling::downscalingFactors[emcdownscaling::kEta]); }); //! cluster pt, assuming m=0 (photons)
} // namespace mincluster

DECLARE_SOA_TABLE(MinClusters, "AOD", "MINCLUSTER", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, skimmedcluster::CollisionId, mincluster::StoredDefinition, mincluster::StoredE, mincluster::StoredEta, mincluster::StoredPhi, mincluster::StoredNCells, mincluster::StoredM02, mincluster::StoredTime,
                  mincluster::Definition<mincluster::StoredDefinition>, mincluster::E<mincluster::StoredE>, mincluster::Eta<mincluster::StoredEta>, mincluster::Phi<mincluster::StoredPhi>, mincluster::NCells<mincluster::StoredNCells>, mincluster::M02<mincluster::StoredM02>, mincluster::Time<mincluster::StoredTime>, mincluster::IsExotic<mincluster::StoredNCells>,
                  mincluster::Pt<mincluster::StoredE, mincluster::StoredEta>);

namespace mintm
{
DECLARE_SOA_INDEX_COLUMN(MinCluster, minCluster); //!
} // namespace mintm

DECLARE_SOA_TABLE(MinMTracks, "AOD", "MINMTRACK", //!
                  mintm::MinClusterId, mincluster::StoredDeltaPhi, mincluster::StoredDeltaEta, mincluster::StoredP, mincluster::StoredPt,
                  mincluster::DeltaPhi<mincluster::StoredDeltaPhi>, mincluster::DeltaEta<mincluster::StoredDeltaEta>, mincluster::TrackP<mincluster::StoredP>, mincluster::TrackPt<mincluster::StoredPt>);

DECLARE_SOA_TABLE(MinMSTracks, "AOD", "MINMSTRACK", //!
                  mintm::MinClusterId, mincluster::StoredDeltaPhi, mincluster::StoredDeltaEta, mincluster::StoredP, mincluster::StoredPt,
                  mincluster::DeltaPhi<mincluster::StoredDeltaPhi>, mincluster::DeltaEta<mincluster::StoredDeltaEta>, mincluster::TrackP<mincluster::StoredP>, mincluster::TrackPt<mincluster::StoredPt>);

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLESREDUX_H_
