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

#ifndef PWGLF_DATAMODEL_REDUCEDLAMBDALAMBDATABLES_H_
#define PWGLF_DATAMODEL_REDUCEDLAMBDALAMBDATABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace redllevent
{
DECLARE_SOA_COLUMN(NumLambda, numLambda, int);     //! Number of lambda
DECLARE_SOA_COLUMN(Centrality, centrality, float); //!
} // namespace redllevent
DECLARE_SOA_TABLE(RedLLEvents, "AOD", "REDLLEVENT",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib,
                  redllevent::Centrality,
                  redllevent::NumLambda);
using RedLLEvent = RedLLEvents::iterator;

namespace lltrack
{
DECLARE_SOA_INDEX_COLUMN(RedLLEvent, redLLEvent);
DECLARE_SOA_COLUMN(LLdId, lldId, int);               //! LL PID
DECLARE_SOA_COLUMN(LLdPx, lldPx, float);             //! LL d Px
DECLARE_SOA_COLUMN(LLdPy, lldPy, float);             //! LL d Py
DECLARE_SOA_COLUMN(LLdPz, lldPz, float);             //! LL d Pz
DECLARE_SOA_COLUMN(LLdx, lldx, float);               //! LL d x
DECLARE_SOA_COLUMN(LLdy, lldy, float);               //! LL d y
DECLARE_SOA_COLUMN(LLdz, lldz, float);               //! LL d z
DECLARE_SOA_COLUMN(LLdMass, lldMass, float);         //! LL d Mass
DECLARE_SOA_COLUMN(LLdd1TPC, lldd1TPC, float);       //! LL dd1 TPC nsigma
DECLARE_SOA_COLUMN(LLdd2TPC, lldd2TPC, float);       //! LL dd2 TPC nsigma
DECLARE_SOA_COLUMN(LLdd1Index, lldd1Index, int64_t); //! LL dd1 global index
DECLARE_SOA_COLUMN(LLdd2Index, lldd2Index, int64_t); //! LL dd2 global index

} // namespace lltrack
DECLARE_SOA_TABLE(LLTracks, "AOD", "LLTRACK",
                  o2::soa::Index<>,
                  lltrack::RedLLEventId,
                  lltrack::LLdId,
                  lltrack::LLdPx,
                  lltrack::LLdPy,
                  lltrack::LLdPz,
                  lltrack::LLdx,
                  lltrack::LLdy,
                  lltrack::LLdz,
                  lltrack::LLdMass,
                  lltrack::LLdd1TPC,
                  lltrack::LLdd2TPC,
                  lltrack::LLdd1Index,
                  lltrack::LLdd2Index);

using LLTrack = LLTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_REDUCEDLAMBDALAMBDATABLES_H_
