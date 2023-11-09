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

// Class for the cluster definition, i.e. what is considered a cluster by the clusterizer.
// The cluster definition contains information about the used algorithm, seed threshold,
// cell energy, gradient as well as timing cut
//
/// \author Florian Jonas <florian.jonas@cern.ch>

#ifndef PWGJE_DATAMODEL_EMCALCLUSTERDEFINITION_H_
#define PWGJE_DATAMODEL_EMCALCLUSTERDEFINITION_H_

#include <string>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
enum class ClusterAlgorithm_t {
  kV1,
  kV3
};
/// \struct EMCALClusterDefinition
/// \brief struct for the cluster definition, i.e. what is considered a cluster by the clusterizer.
/// The cluster definition contains information about the used algorithm, seed threshold,
/// cell energy, gradient as well as timing cut
struct EMCALClusterDefinition {
  ClusterAlgorithm_t algorithm;
  int storageID = -1;              // integer ID used to store cluster definition in a flat table
  int selectedCellType = 1;        // EMCal cell type (CURRENTLY NOT USED TO AVOID MULTIPLE CELL LOOPS)
  std::string name = "kUndefined"; // name of the cluster definition
  double seedEnergy = 0.1;         // seed threshold (GeV)
  double minCellEnergy = 0.05;     // minimum cell energy (GeV)
  double timeMin = -10000;         // minimum time (ns)
  double timeMax = 10000;          // maximum time (ns)
  bool doGradientCut = true;       // apply gradient cut if true
  double gradientCut = -1;         // gradient cut

  // default constructor
  EMCALClusterDefinition() = default;
  // constructor
  EMCALClusterDefinition(ClusterAlgorithm_t pAlgorithm, int pStorageID, int pSelectedCellType, std::string pName, double pSeedEnergy, double pMinCellEnergy, double pTimeMin, double pTimeMax, bool pDoGradientCut, double pGradientCut)
  {
    algorithm = pAlgorithm;
    storageID = pStorageID;
    selectedCellType = pSelectedCellType;
    name = pName;
    seedEnergy = pSeedEnergy;
    minCellEnergy = pMinCellEnergy;
    timeMin = pTimeMin;
    timeMax = pTimeMax;
    doGradientCut = pDoGradientCut;
    gradientCut = pGradientCut;
  }

  // implement comparison operators for int std::string and ClusterAlgorithm_t
  bool operator==(const EMCALClusterDefinition& rhs) const
  {
    return (algorithm == rhs.algorithm && storageID == rhs.storageID && name == rhs.name && seedEnergy == rhs.seedEnergy && minCellEnergy == rhs.minCellEnergy && timeMin == rhs.timeMin && timeMax == rhs.timeMax && gradientCut == rhs.gradientCut && doGradientCut == rhs.doGradientCut);
  }
  bool operator!=(const EMCALClusterDefinition& rhs) const
  {
    return !(*this == rhs);
  }
  bool operator==(const int rhs) const
  {
    return (storageID == rhs);
  }
  bool operator!=(const int rhs) const
  {
    return !(storageID == rhs);
  }
  bool operator==(const std::string& rhs) const
  {
    return (name == rhs);
  }
  bool operator!=(const std::string& rhs) const
  {
    return !(name == rhs);
  }
  bool operator==(const ClusterAlgorithm_t rhs) const
  {
    return (algorithm == rhs);
  }
  bool operator!=(const ClusterAlgorithm_t rhs) const
  {
    return !(algorithm == rhs);
  }

  // cast to int
  operator int() const
  {
    return storageID;
  }

  operator std::string() const
  {
    return name;
  }

  operator ClusterAlgorithm_t() const
  {
    return algorithm;
  }
  std::string toString() const
  {
    return name;
  }
}; // class EMCALClusterDefinition
} // namespace o2::aod
#endif // PWGJE_DATAMODEL_EMCALCLUSTERDEFINITION_H_
