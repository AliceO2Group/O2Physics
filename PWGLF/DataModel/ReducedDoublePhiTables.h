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
/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_
#define PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_

#include <cmath>

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

namespace o2::aod
{
namespace redphievent
{
DECLARE_SOA_COLUMN(NumPos, numPos, int); //! Number of positive Kaon
DECLARE_SOA_COLUMN(NumNeg, numNeg, int); //! Number of negative Kaon
} // namespace redphievent
DECLARE_SOA_TABLE(RedPhiEvents, "AOD", "REDPHIEVENT",
                  o2::soa::Index<>,
                  bc::GlobalBC,
                  bc::RunNumber,
                  timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib,
                  redphievent::NumPos,
                  redphievent::NumNeg);
using RedPhiEvent = RedPhiEvents::iterator;

namespace phitrack
{
DECLARE_SOA_INDEX_COLUMN(RedPhiEvent, redPhiEvent);
DECLARE_SOA_COLUMN(PhiPx, phiPx, float);             //! Phi Px
DECLARE_SOA_COLUMN(PhiPy, phiPy, float);             //! Phi Py
DECLARE_SOA_COLUMN(PhiPz, phiPz, float);             //! Phi Pz
DECLARE_SOA_COLUMN(Phid1Px, phid1Px, float);         //! Phi d1 Px
DECLARE_SOA_COLUMN(Phid1Py, phid1Py, float);         //! Phi d1 Py
DECLARE_SOA_COLUMN(Phid1Pz, phid1Pz, float);         //! Phi d1 Pz
DECLARE_SOA_COLUMN(Phid2Px, phid2Px, float);         //! Phi d2 Px
DECLARE_SOA_COLUMN(Phid2Py, phid2Py, float);         //! Phi d2 Py
DECLARE_SOA_COLUMN(Phid2Pz, phid2Pz, float);         //! Phi d2 Pz
DECLARE_SOA_COLUMN(PhiMass, phiMass, float);         //! Phi Mass
DECLARE_SOA_COLUMN(Phid1Index, phid1Index, int64_t); //! Phi d1 index
DECLARE_SOA_COLUMN(Phid2Index, phid2Index, int64_t); //! Phi d2 index
DECLARE_SOA_COLUMN(Phid1Charge, phid1Charge, float); //! Phi d1 charge
DECLARE_SOA_COLUMN(Phid2Charge, phid2Charge, float); //! Phi d1 charge
DECLARE_SOA_COLUMN(Phid1TPC, phid1TPC, float);       //! TPC nsigma d1
DECLARE_SOA_COLUMN(Phid2TPC, phid2TPC, float);       //! TPC nsigma d2
DECLARE_SOA_COLUMN(Phid1TOFHit, phid1TOFHit, int);   //! TOF hit d1
DECLARE_SOA_COLUMN(Phid2TOFHit, phid2TOFHit, int);   //! TOF hit d2
DECLARE_SOA_COLUMN(Phid1TOF, phid1TOF, float);       //! TOF nsigma d1
DECLARE_SOA_COLUMN(Phid2TOF, phid2TOF, float);       //! TOF nsigma d2
} // namespace phitrack
DECLARE_SOA_TABLE(PhiTracks, "AOD", "PHITRACK",
                  o2::soa::Index<>,
                  phitrack::RedPhiEventId,
                  phitrack::PhiPx,
                  phitrack::PhiPy,
                  phitrack::PhiPz,
                  phitrack::Phid1Px,
                  phitrack::Phid1Py,
                  phitrack::Phid1Pz,
                  phitrack::Phid2Px,
                  phitrack::Phid2Py,
                  phitrack::Phid2Pz,
                  phitrack::PhiMass,
                  phitrack::Phid1Index,
                  phitrack::Phid2Index,
                  phitrack::Phid1Charge,
                  phitrack::Phid2Charge,
                  phitrack::Phid1TPC,
                  phitrack::Phid2TPC,
                  phitrack::Phid1TOFHit,
                  phitrack::Phid2TOFHit,
                  phitrack::Phid1TOF,
                  phitrack::Phid2TOF);

using PhiTrack = PhiTracks::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_REDUCEDDOUBLEPHITABLES_H_
