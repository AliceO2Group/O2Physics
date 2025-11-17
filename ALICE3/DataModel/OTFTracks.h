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
/// \file   OTFTracks.h
/// \author Jesper Karlsson Gumprecht
/// \since  11/11/2025
/// \brief  Table to map track to LUT configuration
///

#ifndef ALICE3_DATAMODEL_OTFTRACKS_H_
#define ALICE3_DATAMODEL_OTFTRACKS_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otftracks
{
DECLARE_SOA_COLUMN(LUTConfigId, lutConfigId, int); //! Index for LUT configuration
} // namespace otftracks

DECLARE_SOA_TABLE(OTFLUTConfigId, "AOD", "OTFLUTConfigId", otftracks::LUTConfigId);
} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFTRACKS_H_
