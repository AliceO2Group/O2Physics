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
/// \file   collisionExtra.h
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table for ALICE 3 collision-related info
///

#ifndef ALICE3_DATAMODEL_PREFILTERDILEPTON_H_
#define ALICE3_DATAMODEL_PREFILTERDILEPTON_H_

// O2 includes
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{

namespace dileptonanalysisflags
{
// DECLARE_SOA_COLUMN(IsMCEventSelected, isMCEventSelected, int);
DECLARE_SOA_COLUMN(IsEventCentSelected, isEventCentSelected, int);
DECLARE_SOA_COLUMN(IsTrackPrefilter, isTrackPrefilter, int);
} // namespace dileptonanalysisflags

// DECLARE_SOA_TABLE(EventMCCuts, "AOD", "EVENTMCCUTS", emanalysisflags::IsMCEventSelected);
DECLARE_SOA_TABLE(DiEventCentCuts, "AOD", "DIEVENTCENTCUTS", dileptonanalysisflags::IsEventCentSelected);
DECLARE_SOA_TABLE(DiTrackPrefilter, "AOD", "DITRACKPREFILTER", dileptonanalysisflags::IsTrackPrefilter);

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_PREFILTERDILEPTON_H_
