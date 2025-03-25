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
namespace fwdtrackrealign
{
// FwdTracksRealign Columns definitions
DECLARE_SOA_INDEX_COLUMN(Collision, collision);    //!
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);      //! FwdTrack index
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t); //! Type of track. See enum ForwardTrackTypeEnum
DECLARE_SOA_COLUMN(X, x, float);                   //! TrackParFwd parameter x
DECLARE_SOA_COLUMN(Y, y, float);                   //! TrackParFwd parameter y
DECLARE_SOA_COLUMN(Z, z, float);                   //! TrackParFwd propagation parameter z
DECLARE_SOA_COLUMN(Phi, phi, float);               //! TrackParFwd parameter phi; (i.e. pt pointing direction)
DECLARE_SOA_COLUMN(Tgl, tgl, float);               //! TrackParFwd parameter tan(\lamba); (\lambda = 90 - \theta_{polar})
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);   //! TrackParFwd parameter: charged inverse transverse momentum; (q/pt)
DECLARE_SOA_COLUMN(IsRemovable, isRemovable, int); //! flag to validate the re-aligned track
DECLARE_SOA_COLUMN(Chi2, chi2, float);             //! Track chi^2

// FwdTracksCovRealign columns definitions
DECLARE_SOA_COLUMN(SigmaX, sigmaX, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaY, sigmaY, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaPhi, sigmaPhi, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaTgl, sigmaTgl, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(RhoXY, rhoXY, int8_t);         //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoPhiX, rhoPhiX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoPhiY, rhoPhiY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglX, rhoTglX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglY, rhoTglY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglPhi, rhoTglPhi, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtX, rho1PtX, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtY, rho1PtY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtPhi, rho1PtPhi, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtTgl, rho1PtTgl, int8_t); //! Covariance matrix in compressed form

// Dynamic and expression columns
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign, //! Sign of the track eletric charge
                           [](float signed1Pt) -> short { return (signed1Pt > 0) ? 1 : -1; });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //!
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //!
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //!
                           [](float pt, float tgl) -> float {
                             return pt * tgl;
                           });

DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //!
                              -1.f * nlog(ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrackrealign::tgl))));
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //!
                              ifnode(nabs(aod::fwdtrackrealign::signed1Pt) < o2::constants::math::Almost0, o2::constants::math::VeryBig, nabs(1.f / aod::fwdtrackrealign::signed1Pt)));
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //!
                              ifnode((nabs(aod::fwdtrackrealign::signed1Pt) < o2::constants::math::Almost0) || (nabs(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrackrealign::tgl)) < o2::constants::math::Almost0), o2::constants::math::VeryBig, 0.5f * (ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrackrealign::tgl)) + 1.f / ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::fwdtrackrealign::tgl))) / nabs(aod::fwdtrackrealign::signed1Pt)));
DECLARE_SOA_EXPRESSION_COLUMN(CXX, cXX, float, //!
                              aod::fwdtrackrealign::sigmaX* aod::fwdtrackrealign::sigmaX);
DECLARE_SOA_EXPRESSION_COLUMN(CXY, cXY, float, //!
                              (aod::fwdtrackrealign::rhoXY / 128.f) * (aod::fwdtrackrealign::sigmaX * aod::fwdtrackrealign::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CYY, cYY, float, //!
                              aod::fwdtrackrealign::sigmaY* aod::fwdtrackrealign::sigmaY);
DECLARE_SOA_EXPRESSION_COLUMN(CPhiX, cPhiX, float, //!
                              (aod::fwdtrackrealign::rhoPhiX / 128.f) * (aod::fwdtrackrealign::sigmaPhi * aod::fwdtrackrealign::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(CPhiY, cPhiY, float, //!
                              (aod::fwdtrackrealign::rhoPhiY / 128.f) * (aod::fwdtrackrealign::sigmaPhi * aod::fwdtrackrealign::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CPhiPhi, cPhiPhi, float, //!
                              aod::fwdtrackrealign::sigmaPhi* aod::fwdtrackrealign::sigmaPhi);
DECLARE_SOA_EXPRESSION_COLUMN(CTglX, cTglX, float, //!
                              (aod::fwdtrackrealign::rhoTglX / 128.f) * (aod::fwdtrackrealign::sigmaTgl * aod::fwdtrackrealign::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(CTglY, cTglY, float, //!
                              (aod::fwdtrackrealign::rhoTglY / 128.f) * (aod::fwdtrackrealign::sigmaTgl * aod::fwdtrackrealign::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CTglPhi, cTglPhi, float, //!
                              (aod::fwdtrackrealign::rhoTglPhi / 128.f) * (aod::fwdtrackrealign::sigmaTgl * aod::fwdtrackrealign::sigmaPhi));
DECLARE_SOA_EXPRESSION_COLUMN(CTglTgl, cTglTgl, float, //!
                              aod::fwdtrackrealign::sigmaTgl* aod::fwdtrackrealign::sigmaTgl);
DECLARE_SOA_EXPRESSION_COLUMN(C1PtY, c1PtY, float, //!
                              (aod::fwdtrackrealign::rho1PtY / 128.f) * (aod::fwdtrackrealign::sigma1Pt * aod::fwdtrackrealign::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtX, c1PtX, float, //!
                              (aod::fwdtrackrealign::rho1PtX / 128.f) * (aod::fwdtrackrealign::sigma1Pt * aod::fwdtrackrealign::sigmaX));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtPhi, c1PtPhi, float, //!
                              (aod::fwdtrackrealign::rho1PtPhi / 128.f) * (aod::fwdtrackrealign::sigma1Pt * aod::fwdtrackrealign::sigmaPhi));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtTgl, c1PtTgl, float, //!
                              (aod::fwdtrackrealign::rho1PtTgl / 128.f) * (aod::fwdtrackrealign::sigma1Pt * aod::fwdtrackrealign::sigmaTgl));
DECLARE_SOA_EXPRESSION_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float, //!
                              aod::fwdtrackrealign::sigma1Pt* aod::fwdtrackrealign::sigma1Pt);
} // namespace fwdtrackrealign

// Tracks including MCH and/or MCH (plus optionally MFT)          //!
DECLARE_SOA_TABLE_FULL(StoredFwdTracksReAlign, "FwdTracksReAlign", "AOD", "FWDTRACKREALIGN",
                       o2::soa::Index<>, fwdtrackrealign::CollisionId, fwdtrackrealign::FwdTrackId, fwdtrackrealign::TrackType, fwdtrackrealign::X, fwdtrackrealign::Y, fwdtrackrealign::Z, fwdtrackrealign::Phi, fwdtrackrealign::Tgl,
                       fwdtrackrealign::Signed1Pt,
                       fwdtrackrealign::Px<fwdtrackrealign::Pt, fwdtrackrealign::Phi>,
                       fwdtrackrealign::Py<fwdtrackrealign::Pt, fwdtrackrealign::Phi>,
                       fwdtrackrealign::Pz<fwdtrackrealign::Pt, fwdtrackrealign::Tgl>,
                       fwdtrackrealign::Sign<fwdtrackrealign::Signed1Pt>,
                       fwdtrackrealign::Chi2,
                       fwdtrackrealign::IsRemovable);

DECLARE_SOA_TABLE_FULL(StoredFwdTrksCovReAlign, "FwdCovsReAlign", "AOD", "FWDCOVREALIGN", //!
                       fwdtrackrealign::SigmaX, fwdtrackrealign::SigmaY, fwdtrackrealign::SigmaPhi, fwdtrackrealign::SigmaTgl, fwdtrackrealign::Sigma1Pt,
                       fwdtrackrealign::RhoXY, fwdtrackrealign::RhoPhiY, fwdtrackrealign::RhoPhiX, fwdtrackrealign::RhoTglX, fwdtrackrealign::RhoTglY,
                       fwdtrackrealign::RhoTglPhi, fwdtrackrealign::Rho1PtX, fwdtrackrealign::Rho1PtY, fwdtrackrealign::Rho1PtPhi, fwdtrackrealign::Rho1PtTgl);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(FwdTracksReAlign, StoredFwdTracksReAlign, "FWDTRKREALIGNEXT", //!
                                fwdtrackrealign::Pt,
                                fwdtrackrealign::Eta,
                                fwdtrackrealign::P); // the table name has here to be the one with EXT which is not nice and under study

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(FwdTrksCovReAlign, StoredFwdTrksCovReAlign, "FWDCOVREALIGNEXT", //!
                                fwdtrackrealign::CXX,
                                fwdtrackrealign::CXY,
                                fwdtrackrealign::CYY,
                                fwdtrackrealign::CPhiX,
                                fwdtrackrealign::CPhiY,
                                fwdtrackrealign::CPhiPhi,
                                fwdtrackrealign::CTglX,
                                fwdtrackrealign::CTglY,
                                fwdtrackrealign::CTglPhi,
                                fwdtrackrealign::CTglTgl,
                                fwdtrackrealign::C1PtX,
                                fwdtrackrealign::C1PtY,
                                fwdtrackrealign::C1PtPhi,
                                fwdtrackrealign::C1PtTgl,
                                fwdtrackrealign::C1Pt21Pt2); // the table name has here to be the one with EXT which is not nice and under study

using FwdTrackRealign = FwdTracksReAlign::iterator;
using FwdTrkCovRealign = FwdTrksCovReAlign::iterator;
using FullFwdTracksRealign = soa::Join<FwdTracksReAlign, FwdTrksCovReAlign>;
using FullFwdTrackRealign = FullFwdTracksRealign::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_FWDTRACKREALIGNTABLES_H_
