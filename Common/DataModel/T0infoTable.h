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
#ifndef O2_T0_INFO_TABLE_H_
#define O2_T0_INFO_TABLE_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace infoFT0
{
DECLARE_SOA_COLUMN(T0A, t0A, float);           //! Collision time A-side, corrected with primary vertex
DECLARE_SOA_COLUMN(T0C, t0C, float);           //! Collision time C-side, corrected with primary vertex
DECLARE_SOA_COLUMN(T0AC, t0AC, float);         //! Collision time (A+C)/2
DECLARE_SOA_COLUMN(T0vertex, t0vertex, float); //! FT0 vertex
} // namespace infoFT0
DECLARE_SOA_TABLE(T0info, "AOD", "T0info",
                  infoFT0::T0A, infoFT0::T0C,
                  infoFT0::T0AC, infoFT0::T0vertex);
using infosFT0 = T0info::iterator;
} // namespace o2::aod

#endif // O2_T0_INFO_TABLE_H_
