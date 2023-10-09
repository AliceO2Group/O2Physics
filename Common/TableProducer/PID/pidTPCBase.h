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
/// \file   pidTPCBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#ifndef COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
#define COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(PIDMults, "AOD", "PIDMults", //! TPC auxiliary table for the PID
                  mult::MultTPC);
using PIDMult = PIDMults::iterator;

} // namespace o2::aod

#endif // COMMON_TABLEPRODUCER_PID_PIDTPCBASE_H_
