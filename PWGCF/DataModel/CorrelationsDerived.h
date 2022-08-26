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
#ifndef O2_ANALYSIS_CFDERIVED_H
#define O2_ANALYSIS_CFDERIVED_H

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
namespace cfcollision
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Centrality/multiplicity value
}
DECLARE_SOA_TABLE_VERSIONED(CFCollisions, "AOD", "CFCOLLISION", 0, //!
                            o2::soa::Index<>,
                            bc::RunNumber, collision::PosZ,
                            cfcollision::Multiplicity, timestamp::Timestamp);
using CFCollision = CFCollisions::iterator;

namespace cftrack
{
DECLARE_SOA_INDEX_COLUMN(CFCollision, cfCollision); //! Index to collision
DECLARE_SOA_COLUMN(Pt, pt, float);                  //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                //! Phi angle
DECLARE_SOA_COLUMN(Sign, sign, int8_t);             //! Sign (positive, negative)
} // namespace cftrack
DECLARE_SOA_TABLE(CFTracks, "AOD", "CFTRACK", //!
                  o2::soa::Index<>,
                  cftrack::CFCollisionId,
                  cftrack::Pt, cftrack::Eta, cftrack::Phi,
                  cftrack::Sign, track::TrackType);
using CFTrack = CFTracks::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_CFDERIVED_H
