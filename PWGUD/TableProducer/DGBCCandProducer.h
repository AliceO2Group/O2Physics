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

#ifndef O2_ANALYSIS_DGBCANDIDATE_H_
#define O2_ANALYSIS_DGBCANDIDATE_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "MathUtils/Utils.h"
#include "Framework/DataTypes.h"
#include "EventFiltering/PWGUD/DGCutparHolder.h"
#include "EventFiltering/PWGUD/DGHelpers.h"
#include "PWGUD/Core/DGMCHelpers.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{

namespace dgbcandidate
{
DECLARE_SOA_COLUMN(BCNum, bcnum, uint64_t);     //! BC number
DECLARE_SOA_COLUMN(FullDiff2, fulldiff2, bool); //! BC contains all tracks if Diff2 event

DECLARE_SOA_INDEX_COLUMN(BC, bc);                      //! BC index
DECLARE_SOA_COLUMN(HasFIT, hasFIT, bool);              //! has FIT signal
DECLARE_SOA_COLUMN(NCollisions, nCollisions, int);     //! number of associated collisions
DECLARE_SOA_COLUMN(NtrwTOF, ntrwTOF, int);             //! number of tracks with TOF
DECLARE_SOA_COLUMN(NPVTracks, nPVTracks, int);         //! number of PV tracks
DECLARE_SOA_COLUMN(NGlobalTracks, nGlobalTracks, int); //! number of global tracks
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, track);          //! track index array

} // namespace dgbcandidate

DECLARE_SOA_TABLE(TrackswTOFInBCs, "AOD", "TRACKSWTOFINBCS", //!
                  o2::soa::Index<>,
                  collision::BCId,
                  dgbcandidate::BCNum,
                  dgbcandidate::TrackIds);

using TrackswTOFInBC = TrackswTOFInBCs::iterator;

DECLARE_SOA_TABLE(IsFullDiff2s, "AOD", "ISFULLDIFF2S", //! To join with TrackswTOFInBCs
                  dgbcandidate::FullDiff2);

using IsFullDiff2 = IsFullDiff2s::iterator;

DECLARE_SOA_TABLE(TracksInBCs, "AOD", "TRACKSINBC", //! Tracks for which a given BC is compatible
                  o2::soa::Index<>,
                  dgbcandidate::BCId,
                  dgbcandidate::HasFIT,
                  dgbcandidate::NCollisions,
                  dgbcandidate::NtrwTOF,
                  dgbcandidate::NPVTracks,
                  dgbcandidate::NGlobalTracks,
                  dgbcandidate::TrackIds);

using TracksInBC = TracksInBCs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_DGBCANDIDATE_H_
