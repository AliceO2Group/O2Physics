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
/// \file   OTFLUT.h
/// \since  23/02/2026
/// \author Jesper & Nicolò
/// \brief  Set of tables for ALICE 3 tracker
///

#ifndef ALICE3_DATAMODEL_OTFLUT_H_
#define ALICE3_DATAMODEL_OTFLUT_H_

#include "ALICE3/Core/DelphesO2TrackSmearer.h"

#include "DataFormatsTOF/CalibLHCphaseTOF.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace otf::lut
{

//   lutHeader_t* mLUTHeader[nLUTs] = {nullptr};
//   lutEntry_t***** mLUTEntry[nLUTs] = {nullptr};

DECLARE_SOA_CCDB_COLUMN(LUTHeader, lutHeader, lutHeader_t, "TOF/Calib/LHCphase"); //!
// DECLARE_SOA_CCDB_COLUMN(LUTHeader, lutHeader, o2::dataformats::CalibLHCphaseTOF, "TOF/Calib/LHCphase"); //!

} // namespace otf::lut

DECLARE_SOA_TIMESTAMPED_TABLE(TOFCalibrationObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "TOFCALIB", //!
                              otf::lut::LUTHeader);
} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFLUT_H_
