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
/// \file   pidTpCBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace pidtpcmult
{
DECLARE_SOA_COLUMN(MultTPC, multTPC, int); //! TPC Multuplicity to evaluate the response resolution
} // namespace pidtpcmult

DECLARE_SOA_TABLE(TPCMult, "AOD", "TPCMult", //! Table of the TPC multiplicity
                  pidtpcmult::MultTPC);
} // namespace o2::aod
