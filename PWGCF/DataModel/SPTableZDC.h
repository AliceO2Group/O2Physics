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

/// \file SPTableZDC.h

#ifndef PWGCF_DATAMODEL_SPTABLEZDC_H_
#define PWGCF_DATAMODEL_SPTABLEZDC_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <cmath>
#include <vector>

namespace o2::aod
{
namespace sptablezdc
{
DECLARE_SOA_COLUMN(Runnumber, runnumber, int);
DECLARE_SOA_COLUMN(Cents, cents, std::vector<float>);
DECLARE_SOA_COLUMN(Vertex, vertex, std::vector<float>);
DECLARE_SOA_COLUMN(Timestamp, timestamp, int64_t);
DECLARE_SOA_COLUMN(QxA, qxA, float);
DECLARE_SOA_COLUMN(QyA, qyA, float);
DECLARE_SOA_COLUMN(QxC, qxC, float);
DECLARE_SOA_COLUMN(QyC, qyC, float);
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);
DECLARE_SOA_COLUMN(EventSelectionFlags, eventSelectionFlags, uint16_t);

} // namespace sptablezdc

DECLARE_SOA_TABLE(SPTableZDC, "AOD", "SPZDC",
                  sptablezdc::Runnumber,
                  sptablezdc::Cents,
                  sptablezdc::Vertex,
                  sptablezdc::Timestamp,
                  sptablezdc::QxA,
                  sptablezdc::QyA,
                  sptablezdc::QxC,
                  sptablezdc::QyC,
                  sptablezdc::IsSelected,
                  sptablezdc::EventSelectionFlags);
} // namespace o2::aod
#endif // PWGCF_DATAMODEL_SPTABLEZDC_H_
