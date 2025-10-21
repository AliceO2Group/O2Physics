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

/// \file FwdTrackReAlignTables.h
/// \brief Table definitions for re-aligned forward tracks
/// \author Chi Zhang <chi.zhang@cern.ch>, CEA-Saclay

#ifndef COMMON_DATAMODEL_FWDTRACKREALIGNTABLES_H_
#define COMMON_DATAMODEL_FWDTRACKREALIGNTABLES_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace fwdtrackrealign
{
DECLARE_SOA_COLUMN(IsRemovable, isRemovable, int); //! flag to check the refit status
}

DECLARE_SOA_TABLE_FULL(StoredFwdTracksReAlign, "FwdTracksReAlign", "AOD", "FWDTRACKREALIGN",
                       o2::soa::Index<>, fwdtrack::CollisionId, fwdtrack::TrackType,
                       fwdtrack::X, fwdtrack::Y, fwdtrack::Z, fwdtrack::Phi, fwdtrack::Tgl,
                       fwdtrack::Signed1Pt, fwdtrack::NClusters, fwdtrack::PDca, fwdtrack::RAtAbsorberEnd,
                       fwdtrackrealign::IsRemovable,
                       fwdtrack::Px<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Py<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Pz<fwdtrack::Pt, fwdtrack::Tgl>,
                       fwdtrack::Sign<fwdtrack::Signed1Pt>,
                       fwdtrack::Chi2, fwdtrack::Chi2MatchMCHMID, fwdtrack::Chi2MatchMCHMFT,
                       fwdtrack::MatchScoreMCHMFT, fwdtrack::MFTTrackId, fwdtrack::MCHTrackId,
                       fwdtrack::MCHBitMap, fwdtrack::MIDBitMap, fwdtrack::MIDBoards,
                       fwdtrack::TrackTime, fwdtrack::TrackTimeRes);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(FwdTracksReAlign, StoredFwdTracksReAlign, "FWDTRKREALIGNEXT", //!
                                fwdtrack::Pt,
                                fwdtrack::Eta,
                                fwdtrack::P);

DECLARE_SOA_TABLE_FULL(StoredFwdTrksCovReAlign, "FwdCovsReAlign", "AOD", "FWDCOVREALIGN",
                       fwdtrack::SigmaX, fwdtrack::SigmaY, fwdtrack::SigmaPhi, fwdtrack::SigmaTgl, fwdtrack::Sigma1Pt,
                       fwdtrack::RhoXY, fwdtrack::RhoPhiY, fwdtrack::RhoPhiX, fwdtrack::RhoTglX, fwdtrack::RhoTglY,
                       fwdtrack::RhoTglPhi, fwdtrack::Rho1PtX, fwdtrack::Rho1PtY, fwdtrack::Rho1PtPhi, fwdtrack::Rho1PtTgl);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(FwdTrksCovReAlign, StoredFwdTrksCovReAlign, "FWDCOVREALIGNEXT", //!
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

using FwdTrackRealign = FwdTracksReAlign::iterator;
using FwdTrkCovRealign = FwdTrksCovReAlign::iterator;
using FullFwdTracksRealign = soa::Join<FwdTracksReAlign, FwdTrksCovReAlign>;
using FullFwdTrackRealign = FullFwdTracksRealign::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_FWDTRACKREALIGNTABLES_H_
