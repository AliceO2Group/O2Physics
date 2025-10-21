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

#ifndef COMMON_DATAMODEL_MFTMCHMATCHINGML_H_
#define COMMON_DATAMODEL_MFTMCHMATCHINGML_H_

#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
DECLARE_SOA_TABLE(FwdTracksML, "AOD", "FWDTRACKML",
                  o2::soa::Index<>,
                  fwdtrack::CollisionId,
                  fwdtrack::TrackType,
                  fwdtrack::X,
                  fwdtrack::Y,
                  fwdtrack::Z,
                  fwdtrack::Phi,
                  fwdtrack::Tgl,
                  fwdtrack::Signed1Pt,
                  fwdtrack::NClusters,
                  fwdtrack::PDca,
                  fwdtrack::RAtAbsorberEnd,
                  fwdtrack::Px<fwdtrack::Pt, fwdtrack::Phi>,
                  fwdtrack::Py<fwdtrack::Pt, fwdtrack::Phi>,
                  fwdtrack::Pz<fwdtrack::Pt, fwdtrack::Tgl>,
                  fwdtrack::Sign<fwdtrack::Signed1Pt>,
                  fwdtrack::Chi2,
                  fwdtrack::Chi2MatchMCHMID,
                  fwdtrack::Chi2MatchMCHMFT,
                  fwdtrack::MatchScoreMCHMFT,
                  fwdtrack::MFTTrackId,
                  fwdtrack::MCHTrackId,
                  fwdtrack::MCHBitMap,
                  fwdtrack::MIDBitMap,
                  fwdtrack::MIDBoards,
                  fwdtrack::TrackTime,
                  fwdtrack::TrackTimeRes,
                  fwdtrack::Eta,
                  fwdtrack::Pt,
                  fwdtrack::P,
                  fwdtrack::FwdDcaX,
                  fwdtrack::FwdDcaY);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_MFTMCHMATCHINGML_H_
