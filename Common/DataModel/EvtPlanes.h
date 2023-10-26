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
DECLARE_SOA_COLUMN(EvtPlUncor, evtPlUncor, float);
DECLARE_SOA_COLUMN(EvtPlRectr, evtPlRectr, float);
DECLARE_SOA_COLUMN(EvtPlTwist, evtPlTwist, float);
DECLARE_SOA_COLUMN(EvtPlFinal, evtPlFinal, float);

DECLARE_SOA_COLUMN(EvtPlBPosUncor, evtPlBPosUncor, float);
DECLARE_SOA_COLUMN(EvtPlBPosRectr, evtPlBPosRectr, float);
DECLARE_SOA_COLUMN(EvtPlBPosTwist, evtPlBPosTwist, float);
DECLARE_SOA_COLUMN(EvtPlBPosFinal, evtPlBPosFinal, float);

DECLARE_SOA_COLUMN(EvtPlBNegUncor, evtPlBNegUncor, float);
DECLARE_SOA_COLUMN(EvtPlBNegRectr, evtPlBNegRectr, float);
DECLARE_SOA_COLUMN(EvtPlBNegTwist, evtPlBNegTwist, float);
DECLARE_SOA_COLUMN(EvtPlBNegFinal, evtPlBNegFinal, float);

DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
} // namespace ep
DECLARE_SOA_TABLE(EvtPlanes, "AOD", "EVTPLANES", //! Table with all event planes.
                  ep::Cent,
                  ep::EvtPlUncor, ep::EvtPlRectr, ep::EvtPlTwist, ep::EvtPlFinal,
                  ep::EvtPlBPosUncor, ep::EvtPlBPosRectr, ep::EvtPlBPosTwist, ep::EvtPlBPosFinal,
                  ep::EvtPlBNegUncor, ep::EvtPlBNegRectr, ep::EvtPlBNegTwist, ep::EvtPlBNegFinal,
                  ep::NTrkBPos, ep::NTrkBNeg);
using EvtPlane = EvtPlanes::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_EVTPLANES_H_
