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

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace fwdtrack
{
// Extra columns for re-aligned forward tracks
DECLARE_SOA_COLUMN(IsRemovable, isRemovable, int); //! flag to validate the re-aligned track
} // namespace fwdtrack

// Tracks including MCH and/or MCH (plus optionally MFT)          //!
DECLARE_SOA_TABLE_FULL(StoredFwdTracksReAlign, "FwdTracksReAlign", "AOD", "FWDTRACKREALIGN",
                       fwdtrack::X, fwdtrack::Y, fwdtrack::Z, fwdtrack::Phi, fwdtrack::Tgl,
                       fwdtrack::Signed1Pt,
                       fwdtrack::Px<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Py<fwdtrack::Pt, fwdtrack::Phi>,
                       fwdtrack::Pz<fwdtrack::Pt, fwdtrack::Tgl>,
                       fwdtrack::Chi2,
                       fwdtrack::IsRemovable);

DECLARE_SOA_EXTENDED_TABLE(FwdTracksReAlign, StoredFwdTracksReAlign, "FWDTRACKREALIGN", //!
                           aod::fwdtrack::Eta,
                           aod::fwdtrack::Pt,
                           aod::fwdtrack::P);

DECLARE_SOA_TABLE_FULL(StoredFwdTrksCovReAlign, "FwdTrksCovReAlign", "AOD", "FWDTRKCOVREALIGN", //!
                       fwdtrack::SigmaX, fwdtrack::SigmaY, fwdtrack::SigmaPhi, fwdtrack::SigmaTgl, fwdtrack::Sigma1Pt,
                       fwdtrack::RhoXY, fwdtrack::RhoPhiY, fwdtrack::RhoPhiX, fwdtrack::RhoTglX, fwdtrack::RhoTglY,
                       fwdtrack::RhoTglPhi, fwdtrack::Rho1PtX, fwdtrack::Rho1PtY, fwdtrack::Rho1PtPhi, fwdtrack::Rho1PtTgl);

DECLARE_SOA_EXTENDED_TABLE(FwdTrksCovReAlign, StoredFwdTrksCovReAlign, "FWDTRKCOVREALIGN", //!
                           aod::fwdtrack::CXX,
                           aod::fwdtrack::CXY,
                           aod::fwdtrack::CYY,
                           aod::fwdtrack::CPhiX,
                           aod::fwdtrack::CPhiY,
                           aod::fwdtrack::CPhiPhi,
                           aod::fwdtrack::CTglX,
                           aod::fwdtrack::CTglY,
                           aod::fwdtrack::CTglPhi,
                           aod::fwdtrack::CTglTgl,
                           aod::fwdtrack::C1PtX,
                           aod::fwdtrack::C1PtY,
                           aod::fwdtrack::C1PtPhi,
                           aod::fwdtrack::C1PtTgl,
                           aod::fwdtrack::C1Pt21Pt2);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_FWDTRACKREALIGNTABLES_H_
