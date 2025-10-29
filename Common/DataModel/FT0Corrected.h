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
#ifndef COMMON_DATAMODEL_FT0CORRECTED_H_
#define COMMON_DATAMODEL_FT0CORRECTED_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace ft0
{
DECLARE_SOA_COLUMN(T0ACorrected, t0ACorrected, float); //! Collision time A-side, corrected with primary vertex
DECLARE_SOA_COLUMN(T0CCorrected, t0CCorrected, float); //! Collision time C-side, corrected with primary vertex
DECLARE_SOA_DYNAMIC_COLUMN(T0AC, t0AC,                 //! Collision time (A+C)/2
                           [](float t0C, float t0A) -> float { return 0.5 * (t0A + t0C); });
DECLARE_SOA_DYNAMIC_COLUMN(T0ACorrectedValid, t0ACorrectedValid, //! Was T0ACorrected computable?
                           [](float t0) -> bool { return t0 < 1e9; });
DECLARE_SOA_DYNAMIC_COLUMN(T0CCorrectedValid, t0CCorrectedValid, //! Was T0CCorrected computable?
                           [](float t0) -> bool { return t0 < 1e9; });
DECLARE_SOA_DYNAMIC_COLUMN(T0resolution, t0resolution, //! Was T0CCorrected computable?
                           [](float t0A, float t0C) -> float { return 0.5f * (t0A - t0C); });
DECLARE_SOA_DYNAMIC_COLUMN(T0ACValid, t0ACValid, //! Was T0AC computable?
                           [](float t0a, float t0c) -> bool { return (t0a < 1e9) && (t0c < 1e9); });

} // namespace ft0
DECLARE_SOA_TABLE(FT0sCorrected, "AOD", "FT0CORRECTED", //! Table with corrected FT0 values
                  ft0::T0ACorrected, ft0::T0CCorrected,
                  ft0::T0AC<ft0::T0ACorrected, ft0::T0CCorrected>,
                  ft0::T0ACorrectedValid<ft0::T0ACorrected>,
                  ft0::T0CCorrectedValid<ft0::T0CCorrected>,
                  ft0::T0ACValid<ft0::T0ACorrected, ft0::T0CCorrected>,
                  ft0::T0resolution<ft0::T0ACorrected, ft0::T0CCorrected>);
using FT0Corrected = FT0sCorrected::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_FT0CORRECTED_H_
