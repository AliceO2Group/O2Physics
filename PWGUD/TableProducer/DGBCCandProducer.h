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
/// \brief
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  30.09.2022

#ifndef PWGUD_TABLEPRODUCER_DGBCCANDPRODUCER_H_
#define PWGUD_TABLEPRODUCER_DGBCCANDPRODUCER_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "MathUtils/Utils.h"

#include <map>
#include <vector>

namespace o2::aod
{

namespace dgbcandidate
{
DECLARE_SOA_COLUMN(BCNum, bcnum, uint64_t); //! BC number

DECLARE_SOA_INDEX_COLUMN(BC, bc);                   //! BC index
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, track);       //! track index array
DECLARE_SOA_ARRAY_INDEX_COLUMN(FwdTrack, fwdtrack); //! fwd track index array

} // namespace dgbcandidate

DECLARE_SOA_TABLE(TracksWGTInBCs, "AOD", "TRKSWGTINBCS", //!
                  o2::soa::Index<>,
                  collision::BCId,
                  bc::RunNumber,
                  dgbcandidate::BCNum,
                  dgbcandidate::TrackIds);

using TracksWGTInBC = TracksWGTInBCs::iterator;

DECLARE_SOA_TABLE(FwdTracksWGTInBCs, "AOD", "FWDTRKSWGTINBCS", //!
                  o2::soa::Index<>,
                  collision::BCId,
                  bc::RunNumber,
                  dgbcandidate::BCNum,
                  dgbcandidate::FwdTrackIds);

using FwdTracksWGTInBC = FwdTracksWGTInBCs::iterator;

} // namespace o2::aod

#endif // PWGUD_TABLEPRODUCER_DGBCCANDPRODUCER_H_
