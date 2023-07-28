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
/// \file   EvtPlanes.h
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Declaration of the table for the (un)corrected Q-vectors for the event plane
/// determination.
///

#ifndef COMMON_DATAMODEL_EVTPLANES_H_
#define COMMON_DATAMODEL_EVTPLANES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace ep
{
DECLARE_SOA_COLUMN(Cent, cent, float);           //! Centrality percentile.
DECLARE_SOA_COLUMN(EvtPlFT0A, evtPlFT0A, float); //! Event plane for FT0A.
DECLARE_SOA_COLUMN(EvtPlFT0C, evtPlFT0C, float); //! Event plane for FT0C.
DECLARE_SOA_COLUMN(EvtPlFV0A, evtPlFV0A, float); //! Event plane for FV0A.
DECLARE_SOA_COLUMN(EvtPlBPos, evtPlBPos, float); //! Event plane for the central barrel, positive pseudorapidity.
DECLARE_SOA_COLUMN(EvtPlBNeg, evtPlBNeg, float); //! Event plane for the central barrel, negative pseudorapidity.

} // namespace ep
DECLARE_SOA_TABLE(EvtPlanes, "AOD", "EVTPLANES", //! Table with all event planes.
                  ep::Cent,
                  ep::EvtPlFT0A, ep::EvtPlFT0C, ep::EvtPlFV0A,
                  ep::EvtPlBPos, ep::EvtPlBNeg);
using EvtPlane = EvtPlanes::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_EVTPLANES_H_
