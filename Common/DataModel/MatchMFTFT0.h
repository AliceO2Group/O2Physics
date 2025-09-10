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

// \file   MatchMFTFT0.h
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief  Declaration of tables useful for the matching of MFT tracks to FT0-C signals
// \date   03/09/24

#ifndef COMMON_DATAMODEL_MATCHMFTFT0_H_
#define COMMON_DATAMODEL_MATCHMFTFT0_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace indices
{                                          // For bctoft0c
DECLARE_SOA_ARRAY_INDEX_COLUMN(FT0, ft0s); // has_ft0s works now, without it doesn't
DECLARE_SOA_ARRAY_INDEX_COLUMN(BC, bcs);   // has_bcs works now, without it doesn't
} // namespace indices
namespace ambii
{ // for MA2T
DECLARE_SOA_INDEX_COLUMN(MFTTrack, track);
} // namespace ambii
DECLARE_SOA_TABLE(MatchedToFT0, "AOD", "MAFT", indices::BCId, indices::FT0Ids);

DECLARE_SOA_TABLE(BCofMFT, "AOD", "BCOFMFT", ambii::MFTTrackId, indices::BCIds);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_MATCHMFTFT0_H_
