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

#ifndef PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_
#define PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(RedCollisions, "AOD", "REDCOLLISION", //! reduced collision table (same structure as the original collision table)
                  o2::soa::Index<>,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::CovXX, collision::CovXY, collision::CovYY, collision::CovXZ, collision::CovYZ, collision::CovZZ,
                  collision::Flags, collision::Chi2, collision::NumContrib,
                  collision::CollisionTime, collision::CollisionTimeRes,
                  bc::RunNumber);

DECLARE_SOA_TABLE(RedPVMults, "AOD", "REDPVMULT", //! Multiplicity from the PV contributors, joinable with reducedCollisions
                  mult::MultNTracksPV);

DECLARE_SOA_TABLE(RedCentFT0Cs, "AOD", "REDCENTFT0C", //! Reduced Run 3 FT0C centrality table, joinable with reducedCollisions
                  cent::CentFT0C);

namespace reducedtracks3body
{
// track parameter definition
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, RedCollisions, ""); //! Collision index
DECLARE_SOA_COLUMN(X, x, float);                                             //!
DECLARE_SOA_COLUMN(Alpha, alpha, float);                                     //!
DECLARE_SOA_COLUMN(Y, y, float);                                             //!
DECLARE_SOA_COLUMN(Z, z, float);                                             //!
DECLARE_SOA_COLUMN(Snp, snp, float);                                         //!
DECLARE_SOA_COLUMN(Tgl, tgl, float);                                         //!
DECLARE_SOA_COLUMN(Signed1Pt, signed1Pt, float);                             //! (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
DECLARE_SOA_EXPRESSION_COLUMN(Phi, phi, float,                               //! Phi of the track, in radians within [0, 2pi)
                              ifnode(nasin(aod::track::snp) + aod::track::alpha < 0.0f, nasin(aod::track::snp) + aod::track::alpha + o2::constants::math::TwoPI,
                                     ifnode(nasin(aod::track::snp) + aod::track::alpha >= o2::constants::math::TwoPI, nasin(aod::track::snp) + aod::track::alpha - o2::constants::math::TwoPI,
                                            nasin(aod::track::snp) + aod::track::alpha)));
DECLARE_SOA_EXPRESSION_COLUMN(Eta, eta, float, //! Pseudorapidity
                              -1.f * nlog(ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl))));
DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! Transverse momentum of the track in GeV/c
                              ifnode(nabs(aod::track::signed1Pt) <= o2::constants::math::Almost0, o2::constants::math::VeryBig, nabs(1.f / aod::track::signed1Pt)));
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign, //! Charge: positive: 1, negative: -1
                           [](float signed1Pt) -> short { return (signed1Pt > 0) ? 1 : -1; });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Momentum in x-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (r * cs - snp * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Momentum in y-direction in GeV/c
                           [](float signed1Pt, float snp, float alpha) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             // FIXME: GCC & clang should optimize to sincosf
                             float cs = cosf(alpha), sn = sinf(alpha);
                             auto r = std::sqrt((1.f - snp) * (1.f + snp));
                             return pt * (snp * cs + r * sn);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Momentum in z-direction in GeV/c
                           [](float signed1Pt, float tgl) -> float {
                             auto pt = 1.f / std::abs(signed1Pt);
                             return pt * tgl;
                           });
DECLARE_SOA_EXPRESSION_COLUMN(P, p, float, //! Momentum in Gev/c
                              ifnode(nabs(aod::track::signed1Pt) <= o2::constants::math::Almost0, o2::constants::math::VeryBig, 0.5f * (ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl)) + 1.f / ntan(o2::constants::math::PIQuarter - 0.5f * natan(aod::track::tgl))) / nabs(aod::track::signed1Pt)));
DECLARE_SOA_DYNAMIC_COLUMN(Rapidity, rapidity, //! Track rapidity, computed under the mass assumption given as input
                           [](float signed1Pt, float tgl, float mass) -> float {
                             const auto pt = 1.f / std::abs(signed1Pt);
                             const auto pz = pt * tgl;
                             const auto p = 0.5f * (std::tan(o2::constants::math::PIQuarter - 0.5f * std::atan(tgl)) + 1.f / std::tan(o2::constants::math::PIQuarter - 0.5f * std::atan(tgl))) * pt;
                             const auto energy = std::sqrt(p * p + mass * mass);
                             return 0.5f * std::log((energy + pz) / (energy - pz));
                           });

// tracks cov matrix parameter definition
DECLARE_SOA_COLUMN(SigmaY, sigmaY, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaZ, sigmaZ, float);        //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaSnp, sigmaSnp, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(SigmaTgl, sigmaTgl, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);    //! Covariance matrix
DECLARE_SOA_COLUMN(RhoZY, rhoZY, int8_t);         //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoSnpY, rhoSnpY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoSnpZ, rhoSnpZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglY, rhoTglY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglZ, rhoTglZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(RhoTglSnp, rhoTglSnp, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtY, rho1PtY, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtZ, rho1PtZ, int8_t);     //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtSnp, rho1PtSnp, int8_t); //! Covariance matrix in compressed form
DECLARE_SOA_COLUMN(Rho1PtTgl, rho1PtTgl, int8_t); //! Covariance matrix in compressed form

DECLARE_SOA_EXPRESSION_COLUMN(CYY, cYY, float, //!
                              aod::track::sigmaY* aod::track::sigmaY);
DECLARE_SOA_EXPRESSION_COLUMN(CZY, cZY, float, //!
                              (aod::track::rhoZY / 128.f) * (aod::track::sigmaZ * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CZZ, cZZ, float, //!
                              aod::track::sigmaZ* aod::track::sigmaZ);
DECLARE_SOA_EXPRESSION_COLUMN(CSnpY, cSnpY, float, //!
                              (aod::track::rhoSnpY / 128.f) * (aod::track::sigmaSnp * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CSnpZ, cSnpZ, float, //!
                              (aod::track::rhoSnpZ / 128.f) * (aod::track::sigmaSnp * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(CSnpSnp, cSnpSnp, float, //!
                              aod::track::sigmaSnp* aod::track::sigmaSnp);
DECLARE_SOA_EXPRESSION_COLUMN(CTglY, cTglY, float, //!
                              (aod::track::rhoTglY / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(CTglZ, cTglZ, float, //!
                              (aod::track::rhoTglZ / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(CTglSnp, cTglSnp, float, //!
                              (aod::track::rhoTglSnp / 128.f) * (aod::track::sigmaTgl * aod::track::sigmaSnp));
DECLARE_SOA_EXPRESSION_COLUMN(CTglTgl, cTglTgl, float, //!
                              aod::track::sigmaTgl* aod::track::sigmaTgl);
DECLARE_SOA_EXPRESSION_COLUMN(C1PtY, c1PtY, float, //!
                              (aod::track::rho1PtY / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaY));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtZ, c1PtZ, float, //!
                              (aod::track::rho1PtZ / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaZ));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtSnp, c1PtSnp, float, //!
                              (aod::track::rho1PtSnp / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaSnp));
DECLARE_SOA_EXPRESSION_COLUMN(C1PtTgl, c1PtTgl, float, //!
                              (aod::track::rho1PtTgl / 128.f) * (aod::track::sigma1Pt * aod::track::sigmaTgl));
DECLARE_SOA_EXPRESSION_COLUMN(C1Pt21Pt2, c1Pt21Pt2, float, //!
                              aod::track::sigma1Pt* aod::track::sigma1Pt);

// tracks extra parameter definition
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, float);                                      //! Momentum at inner wall of the TPC
DECLARE_SOA_COLUMN(Flags, flags, uint32_t);                                                   //! Track flags. Run 2: see TrackFlagsRun2Enum | Run 3: see TrackFlags
DECLARE_SOA_COLUMN(ITSClusterSizes, itsClusterSizes, uint32_t);                               //! Clusters sizes, four bits per a layer, starting from the innermost
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint8_t);                                //! Findable TPC clusters for this track geometry
DECLARE_SOA_COLUMN(TPCNClsFindableMinusFound, tpcNClsFindableMinusFound, int8_t);             //! TPC Clusters: Findable - Found
DECLARE_SOA_COLUMN(TPCNClsFindableMinusCrossedRows, tpcNClsFindableMinusCrossedRows, int8_t); //! TPC Clusters: Findable - crossed rows
DECLARE_SOA_COLUMN(TRDPattern, trdPattern, uint8_t);                                          //! Contributor to the track on TRD layer in bits 0-5, starting from the innermost, bit 6 indicates a potentially split tracklet, bit 7 if the track crossed a padrow
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);                                            //! Chi2 / cluster for the TPC track segment
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, float);                                                  //! Chi2 for the TOF track segment
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                              //! dE/dx signal in the TPC
DECLARE_SOA_COLUMN(TOFExpMom, tofExpMom, float);                                              //! TOF expected momentum obtained in tracking, used to compute the expected times

DECLARE_SOA_EXPRESSION_COLUMN(DetectorMap, detectorMap, uint8_t, //! Detector map version 1, see enum DetectorMapEnum
                              ifnode(aod::track::itsClusterSizes > (uint32_t)0, static_cast<uint8_t>(o2::aod::track::ITS), (uint8_t)0x0) |
                                ifnode(aod::track::tpcNClsFindable > (uint8_t)0, static_cast<uint8_t>(o2::aod::track::TPC), (uint8_t)0x0) |
                                ifnode(aod::track::trdPattern > (uint8_t)0, static_cast<uint8_t>(o2::aod::track::TRD), (uint8_t)0x0) |
                                ifnode((aod::track::tofChi2 >= 0.f) && (aod::track::tofExpMom > 0.f), static_cast<uint8_t>(o2::aod::track::TOF), (uint8_t)0x0));

DECLARE_SOA_DYNAMIC_COLUMN(ITSClsSizeInLayer, itsClsSizeInLayer, //! Size of the ITS cluster in a given layer
                           [](uint32_t itsClusterSizes, int layer) -> uint8_t {
                             if (layer >= 7 || layer < 0) {
                               return 0;
                             }
                             return (itsClusterSizes >> (layer * 4)) & 0xf;
                           });

DECLARE_SOA_DYNAMIC_COLUMN(HasTPC, hasTPC, //! Flag to check if track has a TPC match
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TPC; });
DECLARE_SOA_DYNAMIC_COLUMN(HasTOF, hasTOF, //! Flag to check if track has a TOF measurement
                           [](uint8_t detectorMap) -> bool { return detectorMap & o2::aod::track::TOF; });
DECLARE_SOA_DYNAMIC_COLUMN(IsPVContributor, isPVContributor, //! Run 3: Has this track contributed to the collision vertex fit
                           [](uint8_t flags) -> bool { return (flags & o2::aod::track::PVContributor) == o2::aod::track::PVContributor; });
DECLARE_SOA_DYNAMIC_COLUMN(PIDForTracking, pidForTracking, //! PID hypothesis used during tracking. See the constants in the class PID in PID.h
                           [](uint32_t flags) -> uint32_t { return flags >> 28; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNClsFound, tpcNClsFound, //! Number of found TPC clusters
                           [](uint8_t tpcNClsFindable, int8_t tpcNClsFindableMinusFound) -> int16_t { return (int16_t)tpcNClsFindable - tpcNClsFindableMinusFound; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, //! Number of crossed TPC Rows
                           [](uint8_t tpcNClsFindable, int8_t TPCNClsFindableMinusCrossedRows) -> int16_t { return (int16_t)tpcNClsFindable - TPCNClsFindableMinusCrossedRows; });
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, //! Ratio  crossed rows over findable clusters
                           [](uint8_t tpcNClsFindable, int8_t tpcNClsFindableMinusCrossedRows) -> float {
                             int16_t tpcNClsCrossedRows = (int16_t)tpcNClsFindable - tpcNClsFindableMinusCrossedRows;
                             return (float)tpcNClsCrossedRows / (float)tpcNClsFindable;
                           });

// track PID definition
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the TOF detector for deuteron (recalculated)

} // namespace reducedtracks3body

DECLARE_SOA_TABLE_FULL(StoredRedIUTracks, "RedIUTracks", "AOD", "REDIUTRACK", //! On disk version of the track parameters at inner most update (e.g. ITS) as it comes from the tracking
                       o2::soa::Index<>, reducedtracks3body::CollisionId,
                       reducedtracks3body::X, reducedtracks3body::Alpha,
                       reducedtracks3body::Y, reducedtracks3body::Z, reducedtracks3body::Snp, reducedtracks3body::Tgl,
                       reducedtracks3body::Signed1Pt,
                       // cov matrix
                       reducedtracks3body::SigmaY, reducedtracks3body::SigmaZ, reducedtracks3body::SigmaSnp, reducedtracks3body::SigmaTgl, reducedtracks3body::Sigma1Pt,
                       reducedtracks3body::RhoZY, reducedtracks3body::RhoSnpY, reducedtracks3body::RhoSnpZ, reducedtracks3body::RhoTglY, reducedtracks3body::RhoTglZ,
                       reducedtracks3body::RhoTglSnp, reducedtracks3body::Rho1PtY, reducedtracks3body::Rho1PtZ, reducedtracks3body::Rho1PtSnp, reducedtracks3body::Rho1PtTgl,
                       // tracks extra
                       reducedtracks3body::TPCInnerParam, reducedtracks3body::Flags, reducedtracks3body::ITSClusterSizes,
                       reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusFound, reducedtracks3body::TPCNClsFindableMinusCrossedRows,
                       reducedtracks3body::TRDPattern, reducedtracks3body::TPCChi2NCl, reducedtracks3body::TOFChi2,
                       reducedtracks3body::TPCSignal, reducedtracks3body::TOFExpMom,
                       // TPC PID
                       reducedtracks3body::TPCNSigmaPr, reducedtracks3body::TPCNSigmaPi, reducedtracks3body::TPCNSigmaDe,
                       reducedtracks3body::TOFNSigmaDe,

                       // ----------- dynmaic columns ------------
                       // tracks IU
                       reducedtracks3body::Px<reducedtracks3body::Signed1Pt, reducedtracks3body::Snp, reducedtracks3body::Alpha>,
                       reducedtracks3body::Py<reducedtracks3body::Signed1Pt, reducedtracks3body::Snp, reducedtracks3body::Alpha>,
                       reducedtracks3body::Pz<reducedtracks3body::Signed1Pt, reducedtracks3body::Tgl>,
                       reducedtracks3body::Rapidity<reducedtracks3body::Signed1Pt, reducedtracks3body::Tgl>,
                       reducedtracks3body::Sign<reducedtracks3body::Signed1Pt>,
                       // tracks extra
                       reducedtracks3body::PIDForTracking<reducedtracks3body::Flags>,
                       reducedtracks3body::IsPVContributor<reducedtracks3body::Flags>,
                       reducedtracks3body::HasTPC<reducedtracks3body::DetectorMap>,
                       reducedtracks3body::HasTOF<reducedtracks3body::DetectorMap>,
                       reducedtracks3body::TPCNClsFound<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusFound>,
                       reducedtracks3body::TPCNClsCrossedRows<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusCrossedRows>,
                       reducedtracks3body::ITSClsSizeInLayer<reducedtracks3body::ITSClusterSizes>,
                       reducedtracks3body::TPCCrossedRowsOverFindableCls<reducedtracks3body::TPCNClsFindable, reducedtracks3body::TPCNClsFindableMinusCrossedRows>);

DECLARE_SOA_EXTENDED_TABLE_USER(RedIUTracks, StoredRedIUTracks, "REDIUTRACKEXT", //! Track parameters at inner most update (e.g. ITS) as it comes from the tracking
                                reducedtracks3body::Pt,
                                reducedtracks3body::P,
                                reducedtracks3body::Eta,
                                reducedtracks3body::Phi,
                                // cov matrix
                                reducedtracks3body::CYY,
                                reducedtracks3body::CZY,
                                reducedtracks3body::CZZ,
                                reducedtracks3body::CSnpY,
                                reducedtracks3body::CSnpZ,
                                reducedtracks3body::CSnpSnp,
                                reducedtracks3body::CTglY,
                                reducedtracks3body::CTglZ,
                                reducedtracks3body::CTglSnp,
                                reducedtracks3body::CTglTgl,
                                reducedtracks3body::C1PtY,
                                reducedtracks3body::C1PtZ,
                                reducedtracks3body::C1PtSnp,
                                reducedtracks3body::C1PtTgl,
                                reducedtracks3body::C1Pt21Pt2,
                                // tracks extra
                                reducedtracks3body::DetectorMap);

namespace reduceddecay3body
{
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, RedIUTracks, "_0");       //! Track 0 index
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, RedIUTracks, "_1");       //! Track 1 index
DECLARE_SOA_INDEX_COLUMN_FULL(Track2, track2, int, RedIUTracks, "_2");       //! Track 2 index
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, RedCollisions, ""); //! Collision index
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //! decay3body radius
DECLARE_SOA_COLUMN(Radius, radius, float);                                   //! decay3body phi
DECLARE_SOA_COLUMN(PosZ, posz, float);                                       //! decay3body z position
} // namespace reduceddecay3body

DECLARE_SOA_TABLE(RedDecay3Bodys, "AOD", "REDDECAY3BODY", //! reduced 3-body decay table
                  o2::soa::Index<>, reduceddecay3body::CollisionId, reduceddecay3body::Track0Id, reduceddecay3body::Track1Id, reduceddecay3body::Track2Id);

using ReducedDecay3BodysLinked = soa::Join<RedDecay3Bodys, Decay3BodyDataLink>;
using ReducedDecay3BodyLinked = ReducedDecay3BodysLinked::iterator;

DECLARE_SOA_TABLE(Red3BodyInfo, "AOD", "RED3BODYINFO", //! joinable with RedDecay3Bodys
                  reduceddecay3body::Radius, reduceddecay3body::Phi, reduceddecay3body::PosZ);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_REDUCED3BODYTABLES_H_
