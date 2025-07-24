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
#ifndef PWGLF_DATAMODEL_V0SELECTORTABLES_H_
#define PWGLF_DATAMODEL_V0SELECTORTABLES_H_

#include <cstdint>
#include <Framework/ASoA.h>
namespace o2::aod
{

namespace v0flags
{

enum V0Flags : uint8_t {
  FK0S = 0x1,        // K0S candidate
  FLAMBDA = 0x2,     // Lambda candidate
  FANTILAMBDA = 0x4, // AntiLambda candidate
  FREJECTED = 0x8    // Does not satisfy any of the above, or is randomly rejected
};

DECLARE_SOA_COLUMN(SignalFlag, signalFlag, uint8_t);
DECLARE_SOA_DYNAMIC_COLUMN(IsK0SCandidate, isK0SCandidate, //! Flag to check if V0 is a K0S candidate
                           [](uint8_t flag) -> bool { return flag & o2::aod::v0flags::FK0S; });
DECLARE_SOA_DYNAMIC_COLUMN(IsLambdaCandidate, isLambdaCandidate, //! Flag to check if V0 is a Lambda candidate
                           [](uint8_t flag) -> bool { return flag & o2::aod::v0flags::FLAMBDA; });
DECLARE_SOA_DYNAMIC_COLUMN(IsAntiLambdaCandidate, isAntiLambdaCandidate, //! Flag to check if V0 is a AntiLambda candidate
                           [](uint8_t flag) -> bool { return flag & o2::aod::v0flags::FANTILAMBDA; });
DECLARE_SOA_DYNAMIC_COLUMN(IsRejectedCandidate, isRejectedCandidate, //! Flag to check if V0 is rejected
                           [](uint8_t flag) -> bool { return flag & o2::aod::v0flags::FREJECTED; });
} // namespace v0flags

DECLARE_SOA_TABLE_STAGED(V0SignalFlags, "V0SIGNALFLAGS",
                         v0flags::SignalFlag,
                         v0flags::IsK0SCandidate<v0flags::SignalFlag>,
                         v0flags::IsLambdaCandidate<v0flags::SignalFlag>,
                         v0flags::IsAntiLambdaCandidate<v0flags::SignalFlag>,
                         v0flags::IsRejectedCandidate<v0flags::SignalFlag>);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_V0SELECTORTABLES_H_
