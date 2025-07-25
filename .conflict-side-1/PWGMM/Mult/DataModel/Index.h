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

#ifndef O2_ANALYSIS_INDEX_H_
#define O2_ANALYSIS_INDEX_H_

#include "Framework/AnalysisDataModel.h"
namespace o2::aod
{
namespace idx
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, tracks);
DECLARE_SOA_ARRAY_INDEX_COLUMN(MFTTrack, mfttracks);
} // namespace idx
DECLARE_SOA_TABLE(ParticlesToTracks, "AOD", "P2T", idx::TrackIds);
DECLARE_SOA_TABLE(ParticlesToMftTracks, "AOD", "P2MFTT", idx::MFTTrackIds);
} // namespace o2::aod
#endif // O2_ANALYSIS_INDEX_H_
