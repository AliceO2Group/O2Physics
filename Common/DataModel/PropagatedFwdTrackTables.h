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

/// \file PropagatedFwdTrackTables.h
/// \brief Table definitions for propagated forward tracks
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>
/// \author Luca Micheletti <luca.micheletti@cern.ch>
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef COMMON_DATAMODEL_PROPAGATEDFWDTRACKTABLES_H_
#define COMMON_DATAMODEL_PROPAGATEDFWDTRACKTABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "TrackSelectionTables.h"

namespace o2::aod
{
namespace propfwdtrack
{
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);                                             //! FwdTrack index
DECLARE_SOA_INDEX_COLUMN_FULL(MCHTrack, matchMCHTrack, int, FwdTracks, "_MatchMCHTrack"); //! Index of matched MCH track for GlobalMuonTracks and GlobalForwardTracks
DECLARE_SOA_COLUMN(CXXatDCA, cXXatDCA, float);                                            //! DCAx resolution squared at DCA
DECLARE_SOA_COLUMN(CYYatDCA, cYYatDCA, float);                                            //! DCAy resolution squared at DCA
DECLARE_SOA_COLUMN(CXYatDCA, cXYatDCA, float);                                            //! correlation term of DCAx,y resolution at DCA
DECLARE_SOA_COLUMN(EtaMatchedMCHMID, etaMatchedMCHMID, float);                            //! eta of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(PhiMatchedMCHMID, phiMatchedMCHMID, float);                            //! phi of MCH-MID track in MFT-MCH-MID track at PV
DECLARE_SOA_COLUMN(IsAssociatedToMPC, isAssociatedToMPC, bool);                           //! is assigned to the most probable collision (relevant to TTCA)
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, bool);                                       //! is ambiguous (relevant to TTCA)
} // namespace propfwdtrack

DECLARE_SOA_TABLE_FULL(StoredPropagatedFwdTracks, "PropagatedFwdTracks", "AOD", "PROPFWDTRACK",
                       o2::soa::Index<>, fwdtrack::CollisionId, fwdtrack::TrackType,
                       fwdtrack::X, fwdtrack::Y, fwdtrack::Z, fwdtrack::Phi, fwdtrack::Tgl,
                       fwdtrack::Signed1Pt, fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                       fwdtrack::Px<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Py<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Pz<fwdtrack::Pt, fwdtrack::Tgl>,
                       fwdtrack::Sign<fwdtrack::Signed1Pt>,
                       fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                       fwdtrack::MatchScoreMCHMFT, propfwdtrack::FwdTrackId, fwdtrack::MFTTrackId, propfwdtrack::MCHTrackId,
                       fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                       fwdtrack::TrackTime, fwdtrack::TrackTimeRes, fwdtrack::FwdDcaX, fwdtrack::FwdDcaY,
                       propfwdtrack::CXXatDCA, propfwdtrack::CYYatDCA, propfwdtrack::CXYatDCA,
                       propfwdtrack::EtaMatchedMCHMID, propfwdtrack::PhiMatchedMCHMID,
                       propfwdtrack::IsAssociatedToMPC, propfwdtrack::IsAmbiguous, o2::soa::Marker<1>);

DECLARE_SOA_TABLE_FULL(StoredPropagatedFwdTracksCov, "PropagatedFwdTracksCov", "AOD", "PROPFWDTRACKCOV", //!
                       fwdtrack::SigmaX, fwdtrack::SigmaY, fwdtrack::SigmaPhi, fwdtrack::SigmaTgl, fwdtrack::Sigma1Pt,
                       fwdtrack::RhoXY, fwdtrack::RhoPhiY, fwdtrack::RhoPhiX, fwdtrack::RhoTglX, fwdtrack::RhoTglY,
                       fwdtrack::RhoTglPhi, fwdtrack::Rho1PtX, fwdtrack::Rho1PtY, fwdtrack::Rho1PtPhi, fwdtrack::Rho1PtTgl, o2::soa::Marker<1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(PropagatedFwdTracks, StoredPropagatedFwdTracks, "PROPFWDTRACKEXT", //!
                                fwdtrack::Pt,
                                fwdtrack::Eta,
                                fwdtrack::P);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(PropagatedFwdTracksCov, StoredPropagatedFwdTracksCov, "PROPFWDTRACKCOVEXT", //!
                                fwdtrack::CXX,
                                fwdtrack::CXY,
                                fwdtrack::CYY,
                                fwdtrack::CPhiX,
                                fwdtrack::CPhiY,
                                fwdtrack::CPhiPhi,
                                fwdtrack::CTglX,
                                fwdtrack::CTglY,
                                fwdtrack::CTglPhi,
                                fwdtrack::CTglTgl,
                                fwdtrack::C1PtX,
                                fwdtrack::C1PtY,
                                fwdtrack::C1PtPhi,
                                fwdtrack::C1PtTgl,
                                fwdtrack::C1Pt21Pt2);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PROPAGATEDFWDTRACKTABLES_H_
