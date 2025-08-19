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

/// \file utilsBfieldCCDB.h
/// \brief Utility to set the B field in analysis querying it from CCDB
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova, Italy

#ifndef PWGHF_UTILS_UTILSBFIELDCCDB_H_
#define PWGHF_UTILS_UTILSBFIELDCCDB_H_

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Logger.h>

#include <string> // std::string

/// \brief Sets up the grp object for magnetic field (w/o matCorr for propagation)
/// \param bc is the bunch crossing
/// \param mRunNumber is an int with the run umber of the previous iteration. If at the current iteration it changes, then the grp object is updated
/// \param ccdb is the o2::ccdb::BasicCCDBManager object
/// \param ccdbPathGrp is the path where the GRP oject is stored
/// \param lut is a pointer to the o2::base::MatLayerCylSet object
/// \param isRun2 tells whether we are analysing Run2 converted data or not (different GRP object type)
template <typename TBc>
void initCCDB(TBc const& bc,
              int& mRunNumber,
              o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb,
              std::string const& ccdbPathGrp,
              o2::base::MatLayerCylSet* lut,
              const bool isRun2)
{
  if (mRunNumber != bc.runNumber()) {
    LOGF(info, "====== initCCDB function called (isRun2==%d)", isRun2);
    if (isRun2) { // Run 2 GRP object
      o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
    } else { // Run 3 GRP object
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
    }
    if (lut) {
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    mRunNumber = bc.runNumber();
  }
} /// end initCCDB

#endif // PWGHF_UTILS_UTILSBFIELDCCDB_H_
