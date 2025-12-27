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
/// \author Junlee Kim <junlee.kim@cern.ch>

#ifndef PWGLF_DATAMODEL_REDUCEDHEPTAQUARKTABLES_H_
#define PWGLF_DATAMODEL_REDUCEDHEPTAQUARKTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace redhqevent
{
DECLARE_SOA_COLUMN(NumPhi, numPhi, int);           //! Number of negative K
DECLARE_SOA_COLUMN(NumLambda, numLambda, int);     //! Number of lambda
DECLARE_SOA_COLUMN(Centrality, centrality, float); //!
} // namespace redhqevent
DECLARE_SOA_TABLE(RedHQEvents, "AOD", "REDHQEVENT",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib,
                  redhqevent::Centrality,
                  redhqevent::NumPhi,
                  redhqevent::NumLambda);
using RedHQEvent = RedHQEvents::iterator;

namespace hqtrack
{
DECLARE_SOA_INDEX_COLUMN(RedHQEvent, redHQEvent);
DECLARE_SOA_COLUMN(HQId, hqId, int);               //! HQ ID
DECLARE_SOA_COLUMN(HQPx, hqPx, float);             //! HQ Px
DECLARE_SOA_COLUMN(HQPy, hqPy, float);             //! HQ Py
DECLARE_SOA_COLUMN(HQPz, hqPz, float);             //! HQ Pz
DECLARE_SOA_COLUMN(HQd1Px, hqd1Px, float);         //! HQ d1 Px
DECLARE_SOA_COLUMN(HQd1Py, hqd1Py, float);         //! HQ d1 Py
DECLARE_SOA_COLUMN(HQd1Pz, hqd1Pz, float);         //! HQ d1 Pz
DECLARE_SOA_COLUMN(HQd2Px, hqd2Px, float);         //! HQ d2 Px
DECLARE_SOA_COLUMN(HQd2Py, hqd2Py, float);         //! HQ d2 Py
DECLARE_SOA_COLUMN(HQd2Pz, hqd2Pz, float);         //! HQ d2 Pz
DECLARE_SOA_COLUMN(HQx, hqx, float);               //! HQ x
DECLARE_SOA_COLUMN(HQy, hqy, float);               //! HQ y
DECLARE_SOA_COLUMN(HQz, hqz, float);               //! HQ z
DECLARE_SOA_COLUMN(HQMass, hqMass, float);         //! HQ Mass
DECLARE_SOA_COLUMN(HQd1Index, hqd1Index, int64_t); //! HQ d1 index
DECLARE_SOA_COLUMN(HQd2Index, hqd2Index, int64_t); //! HQ d2 index
DECLARE_SOA_COLUMN(HQd1Charge, hqd1Charge, float); //! HQ d1 charge
DECLARE_SOA_COLUMN(HQd2Charge, hqd2Charge, float); //! HQ d1 charge
DECLARE_SOA_COLUMN(HQd1TPC, hqd1TPC, float);       //! TPC nsigma d1
DECLARE_SOA_COLUMN(HQd2TPC, hqd2TPC, float);       //! TPC nsigma d2
DECLARE_SOA_COLUMN(HQd1TOFHit, hqd1TOFHit, int);   //! TOF hit d1
DECLARE_SOA_COLUMN(HQd2TOFHit, hqd2TOFHit, int);   //! TOF hit d2
DECLARE_SOA_COLUMN(HQd1TOF, hqd1TOF, float);       //! TOF nsigma d1
DECLARE_SOA_COLUMN(HQd2TOF, hqd2TOF, float);       //! TOF nsigma d2

} // namespace hqtrack
DECLARE_SOA_TABLE(HQTracks, "AOD", "HQTRACK",
                  o2::soa::Index<>,
                  hqtrack::RedHQEventId,
                  hqtrack::HQId,
                  hqtrack::HQPx,
                  hqtrack::HQPy,
                  hqtrack::HQPz,
                  hqtrack::HQd1Px,
                  hqtrack::HQd1Py,
                  hqtrack::HQd1Pz,
                  hqtrack::HQd2Px,
                  hqtrack::HQd2Py,
                  hqtrack::HQd2Pz,
                  hqtrack::HQx,
                  hqtrack::HQy,
                  hqtrack::HQz,
                  hqtrack::HQMass,
                  hqtrack::HQd1Index,
                  hqtrack::HQd2Index,
                  hqtrack::HQd1Charge,
                  hqtrack::HQd2Charge,
                  hqtrack::HQd1TPC,
                  hqtrack::HQd2TPC,
                  hqtrack::HQd1TOFHit,
                  hqtrack::HQd2TOFHit,
                  hqtrack::HQd1TOF,
                  hqtrack::HQd2TOF);

using HQTrack = HQTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_REDUCEDHEPTAQUARKTABLES_H_
