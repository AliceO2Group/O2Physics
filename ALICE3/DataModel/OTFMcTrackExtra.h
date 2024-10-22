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
/// \file   OTFMcTrackExtra.h
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>
/// \since  05/08/2024
/// \brief  Table to to hold MC information specifically for xi daughters created in the otf-tracker
///

#ifndef ALICE3_DATAMODEL_OTFMCTRACKEXTRA_H_
#define ALICE3_DATAMODEL_OTFMCTRACKEXTRA_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otf_mctrack_extra
{
DECLARE_SOA_COLUMN(NewPdgCode, newPdgCode, int); //! PDG code (duplicate column but needed for particles created in the otf-tracker)
DECLARE_SOA_COLUMN(IsFromXi, isFromXi, bool);    //! From Xi decayed in otf-tracker
DECLARE_SOA_COLUMN(IsFromL0, isFromL0, bool);    //! From L0 decayed in otf-tracker
} // namespace otf_mctrack_extra
DECLARE_SOA_TABLE(OTFMcExtra, "AOD", "OTFMcExtra",
                  otf_mctrack_extra::NewPdgCode,
                  otf_mctrack_extra::IsFromXi,
                  otf_mctrack_extra::IsFromL0);

using OTFMcTrackExtra = OTFMcExtra::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMCTRACKEXTRA_H_
