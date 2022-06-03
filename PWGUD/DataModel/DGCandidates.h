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
/// \since  03.06.2022

#ifndef O2_ANALYSIS_DGCANDIDATES_H
#define O2_ANALYSIS_DGCANDIDATES_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(DGCandidates, "AOD", "DGCANDIDATES", //!
                  o2::soa::Index<>, bc::RunNumber, timestamp::Timestamp,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib);
using DGCandidate = DGCandidates::iterator;

namespace dgtrack
{
DECLARE_SOA_INDEX_COLUMN(DGCandidate, dgCandidate); //!
DECLARE_SOA_COLUMN(Pt, pt, float);                  //!
DECLARE_SOA_COLUMN(Eta, eta, float);                //!
DECLARE_SOA_COLUMN(Phi, phi, float);                //!
DECLARE_SOA_COLUMN(Sign, sign, int8_t);             //!
} // namespace dgtrack
DECLARE_SOA_TABLE(DGTracks, "AOD", "DGTRACKS", //!
                  o2::soa::Index<>,
                  dgtrack::DGCandidateId,
                  track::Pt, track::Eta, track::Phi,
                  track::Sign);
using DGTrack = DGTracks::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_DGCANDIDATES_H
