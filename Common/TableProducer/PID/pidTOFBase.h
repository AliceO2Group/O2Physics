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
/// \file   pidTOFBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TOF PID tasks.
///

#include "Common/Core/PID/PIDResponse.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace pidtofevtime
{
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);       //! TOF event time
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float); //! TOF event time error
DECLARE_SOA_COLUMN(TOFEvTimeMult, tofEvTimeMult, int); //! TOF event time multiplicity
} // namespace pidtofevtime

DECLARE_SOA_TABLE(TOFEvTime, "AOD", "TOFEvTime", //! Table of the TOF event time
                  pidtofevtime::TOFEvTime,
                  pidtofevtime::TOFEvTimeErr,
                  pidtofevtime::TOFEvTimeMult);
} // namespace o2::aod
