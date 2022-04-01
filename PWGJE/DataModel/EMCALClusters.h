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

// Table definitions for EMCAL analysis clusters
//
// Author: Raymond Ehlers

#ifndef O2_ANALYSIS_DATAMODEL_EMCALCLUSTERS
#define O2_ANALYSIS_DATAMODEL_EMCALCLUSTERS

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace emcalcluster
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                        //!
DECLARE_SOA_COLUMN(ID, id, int);                                       //!
DECLARE_SOA_COLUMN(Energy, energy, float);                             //!
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                     //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //!
DECLARE_SOA_COLUMN(M02, m02, float);                                   //!
DECLARE_SOA_COLUMN(M20, m20, float);                                   //!
DECLARE_SOA_COLUMN(NCells, nCells, int);                               //!
DECLARE_SOA_COLUMN(Time, time, float);                                 //!
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                          //!
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //!
DECLARE_SOA_COLUMN(NLM, nlm, int);                                     //!

} // namespace emcalcluster

DECLARE_SOA_TABLE(EMCALClusters, "AOD", "EMCALCLUSTERS", //!
                  o2::soa::Index<>, emcalcluster::CollisionId, emcalcluster::ID, emcalcluster::Energy,
                  emcalcluster::CoreEnergy, emcalcluster::Eta, emcalcluster::Phi, emcalcluster::M02,
                  emcalcluster::M20, emcalcluster::NCells, emcalcluster::Time,
                  emcalcluster::IsExotic, emcalcluster::DistanceToBadChannel, emcalcluster::NLM);

using EMCALCluster = EMCALClusters::iterator;

} // namespace o2::aod

#endif
