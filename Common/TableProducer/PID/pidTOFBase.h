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

#ifndef COMMON_TABLEPRODUCER_PID_PIDTOFBASE_H_
#define COMMON_TABLEPRODUCER_PID_PIDTOFBASE_H_

#include <string>
#include <vector>

// O2Physics
#include "PID/ParamBase.h"
#include "PID/PIDTOF.h"
#include "Common/DataModel/PIDResponse.h"

static constexpr int nSpecies = 9;
static constexpr int nParameters = 1;
static const std::vector<std::string> particleNames{"El", "Mu", "Pi", "Ka", "Pr", "De", "Tr", "He", "Al"};
static const std::vector<std::string> parameterNames{"Enable"};
static constexpr int defaultParameters[nSpecies][nParameters]{{-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}, {-1}};

namespace o2::aod
{

namespace pidtofevtime
{
// TOF only columns
DECLARE_SOA_COLUMN(UsedForTOFEvTime, usedForTOFEvTime, uint8_t); //! Flag to check if track was used in the TOF event time making
DECLARE_SOA_COLUMN(EvTimeTOF, evTimeTOF, float);                 //! Event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeTOFErr, evTimeTOFErr, float);           //! Error of the event time computed with the TOF detector
DECLARE_SOA_COLUMN(EvTimeTOFMult, evTimeTOFMult, int);           //! Event time multiplicity for TOF
} // namespace pidtofevtime

DECLARE_SOA_TABLE(EvTimeTOFOnly, "AOD", "EvTimeTOFOnly", //! Table for the TOF event time only with TOF. One entry per track.
                  pidtofevtime::UsedForTOFEvTime,
                  pidtofevtime::EvTimeTOF,
                  pidtofevtime::EvTimeTOFErr,
                  pidtofevtime::EvTimeTOFMult);
} // namespace o2::aod

#endif // COMMON_TABLEPRODUCER_PID_PIDTOFBASE_H_
