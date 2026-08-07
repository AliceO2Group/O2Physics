// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoTables.h
/// \brief Datamodel for femto analysis
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_DATAMODEL_FEMTOTABLES_H_
#define PWGCF_FEMTO_DATAMODEL_FEMTOTABLES_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace femtocollisions
{
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::CollisionMaskType);                //! Bitmask for collision selections
DECLARE_SOA_COLUMN(CollisionTag, collisionTag, o2::analysis::femto::datatypes::CollisionTagType); //! Bitmask for collision selections

DECLARE_SOA_COLUMN(PosX, posX, float);                       //! x coordinate of vertex
DECLARE_SOA_COLUMN(PosY, posY, float);                       //! y coordinate of vertex
DECLARE_SOA_COLUMN(PosZ, posZ, float);                       //! z coordinate of vertex
DECLARE_SOA_COLUMN(Mult, mult, float);                       //! Multiplicity estimator set by producer
DECLARE_SOA_COLUMN(Cent, cent, float);                       //! Centrality (~= multiplicity percentile) estimator set by producer
DECLARE_SOA_COLUMN(MagField, magField, int8_t);              //! Magnetic field in kG (5 kG at normal configuration and 2kG in low B field configuration)
DECLARE_SOA_COLUMN(Sphericity, sphericity, float);           //! Sphericity of the event
DECLARE_SOA_COLUMN(Qvec, qvec, float);                       //! qvector
DECLARE_SOA_COLUMN(EventPlaneAngle, eventPlaneAngle, float); //! event plane angle (for corresponding q vector)

namespace lite
{
constexpr float PosZMin = -20.f;
constexpr float PosZMax = 20.f;
constexpr float PosZStep = 0.5f; // bin vtz in 0.5cm steps
constexpr float CentMin = 0.f;
constexpr float CentMax = 100.f;
constexpr float CentStep = 0.5f; // bin centrality in 0.5% steps
constexpr float MultStep = 1.f;  // round multiplicity to nearest integer

inline uint8_t binPosZ(float posZ) { return o2::analysis::femto::utils::binLinear<uint8_t>(posZ, PosZMin, PosZMax, PosZStep); }
inline float unBinPosZ(uint8_t binned) { return o2::analysis::femto::utils::unBinLinear<uint8_t>(binned, PosZMin, PosZStep); }

inline uint8_t binCent(float cent) { return o2::analysis::femto::utils::binLinear<uint8_t>(cent, CentMin, CentMax, CentStep); }
inline float unBinCent(uint8_t binned) { return o2::analysis::femto::utils::unBinLinear<uint8_t>(binned, CentMin, CentStep); }

inline uint16_t binMult(float mult) { return o2::analysis::femto::utils::binLinear<uint16_t>(mult, 0.f, 65535.f, MultStep); } // use full range of uint16_t
inline float unBinMult(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, 0.f, MultStep); }

DECLARE_SOA_COLUMN(BinnedPosZ, binnedPosZ, uint8_t);
DECLARE_SOA_COLUMN(BinnedMult, binnedMult, uint16_t);
DECLARE_SOA_COLUMN(BinnedCent, binnedCent, uint8_t);

DECLARE_SOA_DYNAMIC_COLUMN(PosZ, posZ,
                           [](uint8_t binnedPosZ) -> float {
                             return unBinPosZ(binnedPosZ);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Mult, mult,
                           [](uint16_t binnedMult) -> float {
                             return unBinMult(binnedMult);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Cent, cent,
                           [](uint8_t binnedCent) -> float {
                             return unBinCent(binnedCent);
                           });
} // namespace lite
} // namespace femtocollisions

// table for basic collision information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FCols_001, "FCOL", 1, //! femto collisions
                                   o2::soa::Index<>,
                                   femtocollisions::PosZ,
                                   femtocollisions::Mult,
                                   femtocollisions::Cent,
                                   femtocollisions::MagField);
using FCols = FCols_001;
using FCol = FCols::iterator;
using StoredFCols = StoredFCols_001;

// table for basic collision information, compressed/binned information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteCols_001, "FLITECOL", 1, //! femto collisions, binned information
                                   o2::soa::Index<>,
                                   femtocollisions::lite::BinnedPosZ,
                                   femtocollisions::lite::BinnedMult,
                                   femtocollisions::lite::BinnedCent,
                                   femtocollisions::MagField,
                                   femtocollisions::lite::PosZ<femtocollisions::lite::BinnedPosZ>,
                                   femtocollisions::lite::Mult<femtocollisions::lite::BinnedMult>,
                                   femtocollisions::lite::Cent<femtocollisions::lite::BinnedCent>);
using FLiteCols = FLiteCols_001;
using FLiteCol = FLiteCols::iterator;
using StoredFLiteCols = StoredFLiteCols_001;

// table for collisions selections
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColMasks_001, "FCOLMASK", 1, //! collision masks
                                   femtocollisions::Mask);
using FColMasks = FColMasks_001;
using StoredFColMasks = StoredFColMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FColSphericities_001, "FCOLSPHERICITY", 1, //! sphericity
                                   femtocollisions::Sphericity);
using FColSphericities = FColSphericities_001;

// table for qn values
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColShapes_001, "FCOLSHAPE", 1, //! event shape
                                   femtocollisions::Qvec,
                                   femtocollisions::EventPlaneAngle);
using FColShapes = FColShapes_001;

// table for primary vertex location
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColPos_001, "FCOLPOS", 1, //! full vertex position
                                   femtocollisions::PosX,
                                   femtocollisions::PosY);
using FColPos = FColPos_001;

// table for different multiplicity estimators
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColMults_001, "FCOLMULT", 1,   //! multiplicities
                                   mult::MultFT0A, mult::MultFT0C, //! FIT detectors
                                   mult::MultNTracksPVeta1,        //! number of PV contribs total
                                   mult::MultNTracksPVetaHalf,     //! global track multiplicities
                                   evsel::NumTracksInTimeRange,    //! occupancy (number of track in time range)
                                   evsel::SumAmpFT0CInTimeRange);  //! occupancy (FT0C amplitude in time range)
using FColMults = FColMults_001;

// table for different centrality (multiplicity percentile) estimators
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColCents_001, "FCOLCENT", 1, //! centralities
                                   cent::CentFT0A,               //! centrality from FT0A
                                   cent::CentFT0C,               //! centrality from FT0C
                                   cent::CentFT0M);              //! centrality from FT0M
using FColCents = FColCents_001;

namespace femtobase
// all "basic" information to perform femto analysis, i.e. collision index and kinematics
// split kinematics in stored, i.e. stored in derived data, and dynamic, i.e. can be computed on the fly
{
namespace stored
{
// static columns
DECLARE_SOA_INDEX_COLUMN(FCol, fCol);          //! collision index of femto collision table
DECLARE_SOA_INDEX_COLUMN(FLiteCol, fLiteCol);  //! collision index of femto collision table
DECLARE_SOA_COLUMN(SignedPt, signedPt, float); //! signed pt
DECLARE_SOA_COLUMN(Pt, pt, float);             //! pt
DECLARE_SOA_COLUMN(Eta, eta, float);           //! eta
DECLARE_SOA_COLUMN(Phi, phi, float);           //! phi
DECLARE_SOA_COLUMN(Mass, mass, float);         //! mass of particle
} // namespace stored

namespace dynamic
{
// dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign, //! sign of the track
                           [](float signedPt) -> int {
                             return signedPt > 0.f ? 1 : -1;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! transverse momentum
                           [](float signedPt) -> float {
                             return std::fabs(signedPt);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! momentum in x
                           [](float pt, float phi) -> float {
                             return std::fabs(pt) * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! momentum in y
                           [](float pt, float phi) -> float {
                             return std::fabs(pt) * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! momentum in z
                           [](float pt, float eta) -> float {
                             return std::fabs(pt) * std::sinh(eta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! momentum
                           [](float pt, float eta) -> float {
                             return std::fabs(pt) * std::cosh(eta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, //! theta
                           [](float eta) -> float {
                             return 2.f * std::atan(std::exp(-eta));
                           });
} // namespace dynamic

namespace lite
{

constexpr float EtaMin = -1.5f;
constexpr float EtaMax = 1.5f;
constexpr float EtaStep = (EtaMax - EtaMin) / 65536.f;

constexpr float PhiMin = 0.f;
constexpr float PhiMax = constants::math::TwoPI; // phi is bound on [0,2pi)
constexpr float PhiStep = (PhiMax - PhiMin) / 65536.f;

constexpr float PtMagMin = 0.1f; // lowest pt cut for pions is usually 0.12-0.14 GeV/c
constexpr float PtMagMax = 6.f;  // sensible upper limit for pt

inline uint16_t binEta(float eta) { return o2::analysis::femto::utils::binLinear<uint16_t>(eta, EtaMin, EtaMax, EtaStep); }
inline float unBinEta(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, EtaMin, EtaStep); }

inline uint16_t binPhi(float phi) { return o2::analysis::femto::utils::binLinear<uint16_t>(phi, PhiMin, PhiMax, PhiStep); }
inline float unBinPhi(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, PhiMin, PhiStep); }

inline uint16_t binSignedPt(float signedPt) { return o2::analysis::femto::utils::binLogSigned<uint16_t>(signedPt, PtMagMin, PtMagMax); }
inline float unBinSignedPt(uint16_t binned) { return o2::analysis::femto::utils::unBinLogSigned<uint16_t>(binned, PtMagMin, PtMagMax); }

inline uint16_t binUnsignedPt(float unsignedPt) { return o2::analysis::femto::utils::binLogUnsigned<uint16_t>(unsignedPt, PtMagMin, PtMagMax); }
inline float unBinUnsignedPt(uint16_t binned) { return o2::analysis::femto::utils::unBinLogUnsigned<uint16_t>(binned, PtMagMin, PtMagMax); }

DECLARE_SOA_COLUMN(SignedBinnedPt, signedBinnedPt, uint16_t);
DECLARE_SOA_COLUMN(UnsignedBinnedPt, unsignedBinnedPt, uint16_t);
DECLARE_SOA_COLUMN(BinnedEta, binnedEta, uint16_t);
DECLARE_SOA_COLUMN(BinnedPhi, binnedPhi, uint16_t);

DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,
                           [](uint16_t signedBinnedPt) -> int {
                             return o2::analysis::femto::utils::unBinSign(signedBinnedPt);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(SignedPt, signedPt,
                           [](uint16_t signedBinnedPt) -> float {
                             return unBinSignedPt(signedBinnedPt);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta,
                           [](uint16_t binnedEta) -> float {
                             return unBinEta(binnedEta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi,
                           [](uint16_t binnedPhi) -> float {
                             return unBinPhi(binnedPhi);
                           });

namespace signedpt
{
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](uint16_t signedBinnedPt) -> float {
                             return std::fabs(unBinSignedPt(signedBinnedPt));
                           });
} // namespace signedpt

namespace unsignedpt
{
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt,
                           [](uint16_t unsignedBinnedPt) -> float {
                             return unBinUnsignedPt(unsignedBinnedPt);
                           });
} // namespace unsignedpt

} // namespace lite
} // namespace femtobase

namespace femtotracks
{
// columns for track selections
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::TrackMaskType); //! Bitmask for track selections

// columns for DCA
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                                                                        //! Dca in XY plane
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                                                                          //! Dca in Z direction
DECLARE_SOA_DYNAMIC_COLUMN(Dca, dca, [](float dcaXY, float dcaZ) -> float { return std::hypot(dcaXY, dcaZ); }); //! Dca

// its related information
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);          //! True if track is PV contributer
DECLARE_SOA_COLUMN(ItsNCls, itsNCls, uint8_t);                       //! Number of Its clusters (max 7)
DECLARE_SOA_COLUMN(ItsNClsInnerBarrel, itsNClsInnerBarrel, uint8_t); //! Number of Its clusters in the inner barrel (max 3)
DECLARE_SOA_COLUMN(ItsChi2NCl, itsChi2NCl, float);                   //! Its chi2 / cluster
DECLARE_SOA_COLUMN(ItsClusterSizes, itsClusterSizes, uint32_t);      //! Its cluster sizes (4 bits per layer)

// tpc related information
DECLARE_SOA_COLUMN(TpcSignal, tpcSignal, float);                             //! Tpc signal
DECLARE_SOA_COLUMN(TpcInnerParam, tpcInnerParam, float);                     //! Momentum at inner wall of Tpc
DECLARE_SOA_COLUMN(TpcNClsFound, tpcNClsFound, uint8_t);                     //! Number of Tpc clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, uint8_t);         //! Number of Tpc crossed rows
DECLARE_SOA_DYNAMIC_COLUMN(TpcCrossedRowsOverFound, tpcCrossedRowsOverFound, //! Number of crossed rows over found Tpc clusters
                           [](uint8_t tpcNClsFound, uint8_t tpcNClsCrossedRows) -> float { return static_cast<float>(tpcNClsCrossedRows) / static_cast<float>(tpcNClsFound); });
DECLARE_SOA_COLUMN(TpcNClsShared, tpcNClsShared, uint8_t);         //! Number of shared Tpc clusters
DECLARE_SOA_DYNAMIC_COLUMN(TpcSharedOverFound, tpcSharedOverFound, //! Number of crossed rows over found Tpc clusters
                           [](uint8_t tpcNclsFound, uint8_t tpcNClsShared) -> float { return static_cast<float>(tpcNClsShared) / static_cast<float>(tpcNclsFound); });
DECLARE_SOA_COLUMN(TpcChi2NCl, tpcChi2NCl, float); //! Tpc chi2

// tof related information
DECLARE_SOA_COLUMN(TofBeta, tofBeta, float); //! Tof beta
// tof mass will be stored in mass column

// PID information
// ITS PID information
DECLARE_SOA_COLUMN(ItsNSigmaEl, itsNSigmaEl, float); //! Nsigma separation with the Its for electron
DECLARE_SOA_COLUMN(ItsNSigmaPi, itsNSigmaPi, float); //! Nsigma separation with the Its for pion
DECLARE_SOA_COLUMN(ItsNSigmaKa, itsNSigmaKa, float); //! Nsigma separation with the Its for kaon
DECLARE_SOA_COLUMN(ItsNSigmaPr, itsNSigmaPr, float); //! Nsigma separation with the Its for proton
DECLARE_SOA_COLUMN(ItsNSigmaDe, itsNSigmaDe, float); //! Nsigma separation with the Its for deuteron
DECLARE_SOA_COLUMN(ItsNSigmaTr, itsNSigmaTr, float); //! Nsigma separation with the Its for triton
DECLARE_SOA_COLUMN(ItsNSigmaHe, itsNSigmaHe, float); //! Nsigma separation with the Its for helium3

// TPC PID information
DECLARE_SOA_COLUMN(TpcNSigmaEl, tpcNSigmaEl, float); //! Nsigma separation with the Tpc for electron
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the Tpc for pion
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float); //! Nsigma separation with the Tpc for kaon
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the Tpc for proton
DECLARE_SOA_COLUMN(TpcNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the Tpc for deuteron
DECLARE_SOA_COLUMN(TpcNSigmaTr, tpcNSigmaTr, float); //! Nsigma separation with the Tpc for triton
DECLARE_SOA_COLUMN(TpcNSigmaHe, tpcNSigmaHe, float); //! Nsigma separation with the Tpc for helium3

// TOF PID information
DECLARE_SOA_COLUMN(TofNSigmaEl, tofNSigmaEl, float); //! Nsigma separation with the Tof for electron
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float); //! Nsigma separation with the Tof for pion
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float); //! Nsigma separation with the Tof for kaon
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float); //! Nsigma separation with the Tof for proton
DECLARE_SOA_COLUMN(TofNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the Tof for deuteron
DECLARE_SOA_COLUMN(TofNSigmaTr, tofNSigmaTr, float); //! Nsigma separation with the Tof for triton
DECLARE_SOA_COLUMN(TofNSigmaHe, tofNSigmaHe, float); //! Nsigma separation with the Tof for helium3

DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaEl, tpcitsNSigmaEl, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for electon
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaPi, tpcitsNSigmaPi, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for pion
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaKa, tpcitsNSigmaKa, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for kaon
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaPr, tpcitsNSigmaPr, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for proton
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaDe, tpcitsNSigmaDe, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for deuteron
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaTr, tpcitsNSigmaTr, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for trition
DECLARE_SOA_DYNAMIC_COLUMN(TpcitsNSigmaHe, tpcitsNSigmaHe, [](float tpc, float its) -> float { return std::hypot(tpc, its); }); //! Combined Nsigma separation with Tpc and Its for helium3

DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaEl, tpctofNSigmaEl, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for electons
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaPi, tpctofNSigmaPi, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for pion
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaKa, tpctofNSigmaKa, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for kaon
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaPr, tpctofNSigmaPr, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for proton
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaDe, tpctofNSigmaDe, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for deuteron
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaTr, tpctofNSigmaTr, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for triton
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaHe, tpctofNSigmaHe, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //! Combined Nsigma separation with Tpc and Tof for helium3

} // namespace femtotracks

// table for basic track information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTracks_001, "FTRACK", 1, //! femto tracks
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FTracks = FTracks_001;
using FTrack = FTracks::iterator;
using StoredFTracks = StoredFTracks_001;

// table for basic track information, compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteTracks_001, "FLITETRACK", 1, //! femto tracks, binned kinematics
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>);
using FLiteTracks = FLiteTracks_001;
using FLiteTrack = FLiteTracks::iterator;
using StoredFLiteTracks = StoredFLiteTracks_001;

// table for track selections and PID selections
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTrackMasks_001, "FTRACKMASK", 1, //! track masks
                                   femtotracks::Mask);
using FTrackMasks = FTrackMasks_001;
using StoredFTrackMasks = StoredFTrackMasks_001;

// table for track mass (using tof mass
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTrackMass_001, "FTRACKMASS", 1, //! track mass
                                   femtobase::stored::Mass);
using FTrackMass = FTrackMass_001;
using StoredFTrackMass = StoredFTrackMass_001;

// table for track DCA
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTrackDcas_001, "FTRACKDCAS", 1, //! track dcas
                                   femtotracks::DcaXY,
                                   femtotracks::DcaZ,
                                   femtotracks::Dca<femtotracks::DcaXY, femtotracks::DcaZ>);
using FTrackDcas = FTrackDcas_001;

// table for extra track information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTrackExtras_001, "FTRACKEXTRA", 1, //! track extra information
                                   femtotracks::IsPVContributor,
                                   femtotracks::ItsNCls,
                                   femtotracks::ItsNClsInnerBarrel,
                                   femtotracks::ItsChi2NCl,
                                   femtotracks::ItsClusterSizes,
                                   femtotracks::TpcSignal,
                                   femtotracks::TpcInnerParam,
                                   femtotracks::TpcNClsFound,
                                   femtotracks::TpcNClsCrossedRows,
                                   femtotracks::TpcNClsShared,
                                   femtotracks::TofBeta,
                                   femtotracks::TpcCrossedRowsOverFound<femtotracks::TpcNClsFound, femtotracks::TpcNClsCrossedRows>,
                                   femtotracks::TpcSharedOverFound<femtotracks::TpcNClsFound, femtotracks::TpcNClsShared>);
using FTrackExtras = FTrackExtras_001;

// table for extra PID information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FElectronPids_001, "FELECTRONPID", 1, //! full electron pid
                                   femtotracks::ItsNSigmaEl,
                                   femtotracks::TpcNSigmaEl,
                                   femtotracks::TofNSigmaEl,
                                   femtotracks::TpcitsNSigmaEl<femtotracks::TpcNSigmaEl, femtotracks::ItsNSigmaEl>,
                                   femtotracks::TpctofNSigmaEl<femtotracks::TpcNSigmaEl, femtotracks::TofNSigmaEl>);
using FElectronPids = FElectronPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FPionPids_001, "FPIONPID", 1, //! full pion pid
                                   femtotracks::ItsNSigmaPi,
                                   femtotracks::TpcNSigmaPi,
                                   femtotracks::TofNSigmaPi,
                                   femtotracks::TpcitsNSigmaPi<femtotracks::TpcNSigmaPi, femtotracks::ItsNSigmaPi>,
                                   femtotracks::TpctofNSigmaPi<femtotracks::TpcNSigmaPi, femtotracks::TofNSigmaPi>);
using FPionPids = FPionPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FKaonPids_001, "FKAONPID", 1, //! full kaon pid
                                   femtotracks::ItsNSigmaKa,
                                   femtotracks::TpcNSigmaKa,
                                   femtotracks::TofNSigmaKa,
                                   femtotracks::TpcitsNSigmaKa<femtotracks::TpcNSigmaKa, femtotracks::ItsNSigmaKa>,
                                   femtotracks::TpctofNSigmaKa<femtotracks::TpcNSigmaKa, femtotracks::TofNSigmaKa>);
using FKaonPids = FKaonPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FProtonPids_001, "FPROTONPID", 1, //! full proton pid
                                   femtotracks::ItsNSigmaPr,
                                   femtotracks::TpcNSigmaPr,
                                   femtotracks::TofNSigmaPr,
                                   femtotracks::TpcitsNSigmaPr<femtotracks::TpcNSigmaPr, femtotracks::ItsNSigmaPr>,
                                   femtotracks::TpctofNSigmaPr<femtotracks::TpcNSigmaPr, femtotracks::TofNSigmaPr>);
using FProtonPids = FProtonPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FDeuteronPids_001, "FDEUTERONPID", 1, //! full deuteron pid
                                   femtotracks::ItsNSigmaDe,
                                   femtotracks::TpcNSigmaDe,
                                   femtotracks::TofNSigmaDe,
                                   femtotracks::TpcitsNSigmaDe<femtotracks::TpcNSigmaDe, femtotracks::ItsNSigmaDe>,
                                   femtotracks::TpctofNSigmaDe<femtotracks::TpcNSigmaDe, femtotracks::TofNSigmaDe>);
using FDeuteronPids = FDeuteronPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTritonPids_001, "FTRITONPID", 1, //! full triton pid
                                   femtotracks::ItsNSigmaTr,
                                   femtotracks::TpcNSigmaTr,
                                   femtotracks::TofNSigmaTr,
                                   femtotracks::TpcitsNSigmaTr<femtotracks::TpcNSigmaTr, femtotracks::ItsNSigmaTr>,
                                   femtotracks::TpctofNSigmaTr<femtotracks::TpcNSigmaTr, femtotracks::TofNSigmaTr>);
using FTritonPids = FTritonPids_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FHeliumPids_001, "FHELIUMPID", 1, //! full helium3 pid
                                   femtotracks::ItsNSigmaHe,
                                   femtotracks::TpcNSigmaHe,
                                   femtotracks::TofNSigmaHe,
                                   femtotracks::TpcitsNSigmaHe<femtotracks::TpcNSigmaHe, femtotracks::ItsNSigmaHe>,
                                   femtotracks::TpctofNSigmaHe<femtotracks::TpcNSigmaHe, femtotracks::TofNSigmaHe>);
using FHeliumPids = FHeliumPids_001;

using FTrackPids = soa::Join<FElectronPids, FPionPids, FKaonPids, FProtonPids, FDeuteronPids, FTritonPids, FHeliumPids>;

namespace femtotwotrackresonances
{
// columns for resonance bit masks
// DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::TwoTrackResonanceMaskType); //! Bitmask for resonance selections
DECLARE_SOA_COLUMN(MaskPosDau, maskPosDau, o2::analysis::femto::datatypes::TrackMaskType); //! Bitmask for positive daughter
DECLARE_SOA_COLUMN(PosDauHasHighMomentum, posDauHasHighMomentum, bool);                    //! switch for pid threshold
DECLARE_SOA_COLUMN(MaskNegDau, maskNegDau, o2::analysis::femto::datatypes::TrackMaskType); //! Bitmask for negative daughter
DECLARE_SOA_COLUMN(NegDauHasHighMomentum, negDauHasHighMomentum, bool);                    //! switch for pid threshold

// id columns for resonance daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FTracks, "_PosDau"); //! index column for positive daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FTracks, "_NegDau"); //! index column for negative daughter track
} // namespace femtotwotrackresonances
// table for phis
DECLARE_SOA_TABLE_STAGED_VERSIONED(FPhis_001, "FPHI", 1, //! femto phis
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauId,
                                   femtotwotrackresonances::NegDauId,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FPhis = FPhis_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FPhiMasks_001, "FPHIMASK", 1, //! mask for phis
                                   femtotwotrackresonances::MaskPosDau,
                                   femtotwotrackresonances::PosDauHasHighMomentum,
                                   femtotwotrackresonances::MaskNegDau,
                                   femtotwotrackresonances::NegDauHasHighMomentum);
using FPhiMasks = FPhiMasks_001;

// table for kstars
DECLARE_SOA_TABLE_STAGED_VERSIONED(FKstar0s_001, "FKSTAR0", 1, //! femto k0star
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, //! +1 for k0star and -1 for k0starbar
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauId,
                                   femtotwotrackresonances::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FKstar0s = FKstar0s_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FKstar0Masks_001, "FKSTAR0MASK", 1, //! k0star masks
                                   femtotwotrackresonances::MaskPosDau,
                                   femtotwotrackresonances::PosDauHasHighMomentum,
                                   femtotwotrackresonances::MaskNegDau,
                                   femtotwotrackresonances::NegDauHasHighMomentum);
using FKstar0Masks = FKstar0Masks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FRho0s_001, "FRHO0", 1, //! femto rho0s
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauId,
                                   femtotwotrackresonances::NegDauId,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FRho0s = FRho0s_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FRho0Masks_001, "FRHO0MASK", 1, //! rho0s masks
                                   femtotwotrackresonances::MaskPosDau,
                                   femtotwotrackresonances::PosDauHasHighMomentum,
                                   femtotwotrackresonances::MaskNegDau,
                                   femtotwotrackresonances::NegDauHasHighMomentum);
using FRho0Masks = FRho0Masks_001;

namespace femtov0s
{
// columns for bit masks
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::V0MaskType); //! Bitmask for v0 selections
//
namespace legacy001
{
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::V0MaskType001); //! Bitmask for v0 selections
}

// columns for debug information
DECLARE_SOA_COLUMN(MassAnti, massAnti, float);                 //! mass of particle using antiparticle hypothesis (for Lambda/AntiLambda extra table)
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);             //! Mass of Lambda (for k0short table)
DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float);     //! Mass of AntiLambda (for k0short table)
DECLARE_SOA_COLUMN(MassK0short, massK0short, float);           //! Mass of K0short (for lambda/antitlambda table)
DECLARE_SOA_COLUMN(CosPa, cosPa, float);                       //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(DauDca, dauDca, float);                     //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(StrangeTofPosDau, strangeTofPosDau, float); //! TOF Strangeness for positive daughter
DECLARE_SOA_COLUMN(StrangeTofNegDau, strangeTofNegDau, float); //! TOF Strangeness for negative daughter
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);           //! Lambda transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);               //! x coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);               //! y coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);               //! z coordinate of Lambda decay vertex
DECLARE_SOA_DYNAMIC_COLUMN(DecayVtx, decayVtx,                 //! distance of decay vertex from nominal interaction point
                           [](float vtxX, float vtxY, float vtxZ) -> float {
                             return std::hypot(vtxX, vtxY, vtxZ);
                           });

// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FTracks, "_PosDau"); //! index column for positive daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FTracks, "_NegDau"); //! index column for negative daughter track

namespace lite
{

constexpr float LambdaMassMin = 1.05f; // Kinematic lower limit
constexpr float LambdaMassMax = 1.30f;
constexpr float LambdaMassStep = (LambdaMassMax - LambdaMassMin) / 65536.f;

// K0short: PDG mass 0.497611 GeV roughly ±100 MeV window
constexpr float K0shortMassMin = 0.4f;
constexpr float K0shortMassMax = 0.6f;
constexpr float K0shortMassStep = (K0shortMassMax - K0shortMassMin) / 65536.f;

inline uint16_t binLambdaMass(float lambdaMass) { return o2::analysis::femto::utils::binLinear<uint16_t>(lambdaMass, LambdaMassMin, LambdaMassMax, LambdaMassStep); }
inline float unBinLambdaMass(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, LambdaMassMin, LambdaMassStep); }

inline uint16_t binK0shortMass(float lambdaMass) { return o2::analysis::femto::utils::binLinear<uint16_t>(lambdaMass, K0shortMassMin, K0shortMassMax, K0shortMassStep); }
inline float unBinK0shortMass(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, K0shortMassMin, K0shortMassStep); }

DECLARE_SOA_COLUMN(BinnedLambdaMass, binnedLambdaMass, uint16_t);
DECLARE_SOA_COLUMN(BinnedK0shortMass, binnedK0shortMass, uint16_t);

DECLARE_SOA_DYNAMIC_COLUMN(LambdaMass, lambdaMass,
                           [](uint16_t binnedLambdaMass) -> float {
                             return unBinLambdaMass(binnedLambdaMass);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(K0shortMass, k0shortMass,
                           [](uint16_t binnedK0shortMass) -> float {
                             return unBinK0shortMass(binnedK0shortMass);
                           });

// index for lite track table
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FLiteTracks, "_PosDau"); //! index column for positive daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FLiteTracks, "_NegDau"); //! index column for negative daughter track

} // namespace lite
} // namespace femtov0s

// table for basic lambda information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdas_001, "FLAMBDA", 1, //! femto lambdas
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, // use sign to differentiate between lambda (+1) and antilambda (-1)
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass, // mass of the lambda/antilambda depending on the sign of the pt
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FLambdas = FLambdas_001;
using StoredFLambdas = StoredFLambdas_001;

// table for basic lambda information with compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteLambdas_001, "FLITELAMBDA", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtov0s::lite::BinnedLambdaMass,
                                   femtov0s::lite::PosDauId,
                                   femtov0s::lite::NegDauId,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtov0s::lite::LambdaMass<femtov0s::lite::BinnedLambdaMass>);
using FLiteLambdas = FLiteLambdas_001;
using FLiteLambda = FLiteLambdas::iterator;
using StoredFLiteLambdas = StoredFLiteLambdas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdaMasks_001, "FLAMBDAMASK", 1, //! legacy lambda mask
                                   femtov0s::legacy001::Mask);

DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdaMasks_002, "FLAMBDAMASK", 2, //! lambda mask
                                   femtov0s::Mask);

using FLambdaMasks = FLambdaMasks_002;
using StoredFLambdaMasks = StoredFLambdaMasks_002;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdaExtras_001, "FLAMBDAEXTRA", 1, //! lambda extra information
                                   femtov0s::MassAnti,                   // put mass of antiparticle, i.e. antilambda mass for lambdas and vice versa
                                   femtov0s::MassK0short,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
                                   femtov0s::StrangeTofPosDau,
                                   femtov0s::StrangeTofNegDau,
                                   femtov0s::TransRadius,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::DecayVtx<femtov0s::DecayVtxX, femtov0s::DecayVtxY, femtov0s::DecayVtxZ>);

using FLambdaExtras = FLambdaExtras_001;

// table for basic k0short information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shorts_001, "FK0SHORT", 1, //! femto k0shorts
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FK0shorts = FK0shorts_001;
using StoredFK0shorts = StoredFK0shorts_001;

// table for basic lambda information with compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteK0shorts_001, "FLITEK0SHORT", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::UnsignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtov0s::lite::BinnedK0shortMass,
                                   femtov0s::lite::PosDauId,
                                   femtov0s::lite::NegDauId,
                                   femtobase::lite::unsignedpt::Pt<femtobase::lite::UnsignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtov0s::lite::K0shortMass<femtov0s::lite::BinnedK0shortMass>);
using FLiteK0shorts = FLiteK0shorts_001;
using FLiteK0short = FLiteK0shorts::iterator;
using StoredFLiteK0shorts = StoredFLiteK0shorts_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shortMasks_001, "FK0SHORTMASK", 1, //! legacy k0short masks
                                   femtov0s::legacy001::Mask);

DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shortMasks_002, "FK0SHORTMASK", 2, //! k0short masks
                                   femtov0s::Mask);
using FK0shortMasks = FK0shortMasks_002;
using StoredFK0shortMasks = StoredFK0shortMasks_002;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shortExtras_001, "FK0SHORTEXTRA", 1, //! k0short extra information
                                   femtov0s::MassLambda,
                                   femtov0s::MassAntiLambda,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
                                   femtov0s::StrangeTofPosDau,
                                   femtov0s::StrangeTofNegDau,
                                   femtov0s::TransRadius,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::DecayVtx<femtov0s::DecayVtxX, femtov0s::DecayVtxY, femtov0s::DecayVtxZ>);

using FK0shortExtras = FK0shortExtras_001;

namespace femtokinks
{
// columns for bit masks
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::KinkMaskType); //! Bitmask for kink selections

// columns for debug information
DECLARE_SOA_COLUMN(KinkAngle, kinkAngle, float);     //! Kink angle between mother and charged daughter at decay vertex
DECLARE_SOA_COLUMN(DcaMothToPV, dcaMothToPV, float); //! DCA of the mother track to the primary vertex
DECLARE_SOA_COLUMN(DcaDaugToPV, dcaDaugToPV, float); //! DCA of the charged daughter track to the primary vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of decay vertex (relative to PV)
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of decay vertex (relative to PV)
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of decay vertex (relative to PV)
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Transverse decay radius from PV

// id column for charged daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(ChaDau, chaDau, int32_t, FTracks, "_ChaDau"); //!

namespace lite
{
// Sigma (using a shared window covering both Sigma- 1.19745 and Sigma+ 1.18937) with roughly +-100 MeV
constexpr float SigmaMassMin = 1.1f;
constexpr float SigmaMassMax = 1.3f;
constexpr float SigmaMassStep = (SigmaMassMax - SigmaMassMin) / 65536.f;

inline uint16_t binSigmaMass(float mass) { return o2::analysis::femto::utils::binLinear<uint16_t>(mass, SigmaMassMin, SigmaMassMax, SigmaMassStep); }
inline float unBinSigmaMass(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, SigmaMassMin, SigmaMassStep); }

DECLARE_SOA_COLUMN(BinnedSigmaMass, binnedSigmaMass, uint16_t);

DECLARE_SOA_DYNAMIC_COLUMN(SigmaMass, sigmaMass,
                           [](uint16_t binnedSigmaMass) -> float {
                             return unBinSigmaMass(binnedSigmaMass);
                           });

// index for lite track table
DECLARE_SOA_INDEX_COLUMN_FULL(ChaDau, chaDau, int32_t, FLiteTracks, "_ChaDau"); //! charged daughter index, into FLiteTracks
} // namespace lite

} // namespace femtokinks

// table for basic sigma information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmas_001, "FSIGMA", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtokinks::ChaDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);

// table for basic sigma information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmas_002, "FSIGMA", 2,
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtokinks::ChaDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FSigmas = FSigmas_002;
using FSigma = FSigmas::iterator;
using StoredFSigmas = StoredFSigmas_002;

// table for basic sigma information, compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteSigmas_001, "FLITESIGMA", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtokinks::lite::BinnedSigmaMass,
                                   femtokinks::lite::ChaDauId,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtokinks::lite::SigmaMass<femtokinks::lite::BinnedSigmaMass>);
using FLiteSigmas = FLiteSigmas_001;
using FLiteSigma = FLiteSigmas::iterator;
using StoredFLiteSigmas = StoredFLiteSigmas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaMasks_001, "FSIGMAMASKS", 1,
                                   femtokinks::Mask);
using FSigmaMasks = FSigmaMasks_001;
using FSigmamask = FSigmaMasks::iterator;
using StoredFSigmaMasks = StoredFSigmaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaExtras_001, "FSIGMAEXTRAS", 1,
                                   femtokinks::KinkAngle,
                                   femtokinks::DcaDaugToPV,
                                   femtokinks::DcaMothToPV,
                                   femtokinks::DecayVtxX,
                                   femtokinks::DecayVtxY,
                                   femtokinks::DecayVtxZ,
                                   femtokinks::TransRadius);

using FSigmaExtras = FSigmaExtras_001;

// table for basic sigma plus information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaPlus_001, "FSIGMAPLUS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtokinks::ChaDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FSigmaPlus = FSigmaPlus_001;
using FSigmaPlusRow = FSigmaPlus::iterator;
using StoredFSigmaPlus = StoredFSigmaPlus_001;

// table for basic sigma plus information, compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteSigmaPlus_001, "FLITESIGMAPLUS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtokinks::lite::BinnedSigmaMass,
                                   femtokinks::lite::ChaDauId,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtokinks::lite::SigmaMass<femtokinks::lite::BinnedSigmaMass>);
using FLiteSigmaPlus = FLiteSigmaPlus_001;
using FLiteSigmaPlusRow = FLiteSigmaPlus::iterator;
using StoredFLiteSigmaPlus = StoredFLiteSigmaPlus_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaPlusMasks_001, "FSIGMAPLUSMASKS", 1,
                                   femtokinks::Mask);
using FSigmaPlusMasks = FSigmaPlusMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaPlusExtras_001, "FSIGMAPLUSEXTRAS", 1,
                                   femtokinks::KinkAngle,
                                   femtokinks::DcaDaugToPV,
                                   femtokinks::DcaMothToPV,
                                   femtokinks::DecayVtxX,
                                   femtokinks::DecayVtxY,
                                   femtokinks::DecayVtxZ,
                                   femtokinks::TransRadius);

using FSigmaPlusExtras = FSigmaPlusExtras_001;

namespace femtocascades
{
// columns for cascade bit masks
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::CascadeMaskType); //! Bitmask for cascade selections

namespace legacy001
{
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::CascadeMaskType001); //! Bitmask for cascade selections
}

// columns for cascad debug information
DECLARE_SOA_COLUMN(MassXi, massXi, float);                         //! Mass of xi
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);                   //! Mass of omega
DECLARE_SOA_COLUMN(CascadeCosPa, cascadeCosPa, float);             //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(CascadeDauDca, cascadeDauDca, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(CascadeTransRadius, cascadeTransRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);                 //! Lambda daughter mass
DECLARE_SOA_COLUMN(LambdaCosPa, lambdaCosPa, float);               //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(LambdaDauDca, lambdaDauDca, float);             //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(LambdaTransRadius, lambdaTransRadius, float);   //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaDcaToPv, lambdaDcaToPv, float);           //! Lambda transvers radius
DECLARE_SOA_COLUMN(StrangeTofBachelor, strangeTofBachelor, float); //! Lambda transvers radius

// id columns for bachelor
// following same style as strangeness tables were we do not store the id of the lambda, but its daughters
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int32_t, FTracks, "_Bachelor"); //! bachelor id

namespace lite
{

// Xi-: PDG mass 1.32171 GeV, roughly ±100 MeV window
constexpr float XiMassMin = 1.22f;
constexpr float XiMassMax = 1.42f;
constexpr float XiMassStep = (XiMassMax - XiMassMin) / 65536.f;

// Omega-: PDG mass 1.67245 GeV, roughly ±100 MeV window
constexpr float OmegaMassMin = 1.57f;
constexpr float OmegaMassMax = 1.77f;
constexpr float OmegaMassStep = (OmegaMassMax - OmegaMassMin) / 65536.f;

inline uint16_t binXiMass(float mass) { return o2::analysis::femto::utils::binLinear<uint16_t>(mass, XiMassMin, XiMassMax, XiMassStep); }
inline float unBinXiMass(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, XiMassMin, XiMassStep); }

inline uint16_t binOmegaMass(float mass) { return o2::analysis::femto::utils::binLinear<uint16_t>(mass, OmegaMassMin, OmegaMassMax, OmegaMassStep); }
inline float unBinOmegaMass(uint16_t binned) { return o2::analysis::femto::utils::unBinLinear<uint16_t>(binned, OmegaMassMin, OmegaMassStep); }

DECLARE_SOA_COLUMN(BinnedXiMass, binnedXiMass, uint16_t);
DECLARE_SOA_COLUMN(BinnedOmegaMass, binnedOmegaMass, uint16_t);

DECLARE_SOA_DYNAMIC_COLUMN(XiMass, xiMass,
                           [](uint16_t binnedXiMass) -> float {
                             return unBinXiMass(binnedXiMass);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(OmegaMass, omegaMass,
                           [](uint16_t binnedOmegaMass) -> float {
                             return unBinOmegaMass(binnedOmegaMass);
                           });

// index for lite track table
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int32_t, FLiteTracks, "_Bachelor"); //! bachelor index, into FLiteTracks
} // namespace lite
} // namespace femtocascades

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXis_001, "FXI", 1, //! femto xis
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtocascades::BachelorId,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FXis = FXis_001;
using FXi = FXis::iterator;
using StoredFXis = StoredFXis_001;

// table for basic xi information, compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteXis_001, "FLITEXI", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtocascades::lite::BinnedXiMass,
                                   femtocascades::lite::BachelorId,
                                   femtov0s::lite::PosDauId,
                                   femtov0s::lite::NegDauId,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtocascades::lite::XiMass<femtocascades::lite::BinnedXiMass>);
using FLiteXis = FLiteXis_001;
using FLiteXi = FLiteXis::iterator;
using StoredFLiteXis = StoredFLiteXis_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXiMasks_001, "FXIMASK", 1, //! legacy xi masks
                                   femtocascades::legacy001::Mask);

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXiMasks_002, "FXIMASK", 2, //! xi masks
                                   femtocascades::Mask);
using FXiMasks = FXiMasks_002;
using StoredFXiMasks = StoredFXiMasks_002;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXiExtras_001, "FXIEXTRA", 1, //! xi extra information
                                   femtocascades::MassOmega,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaMass,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv,
                                   femtocascades::StrangeTofBachelor,
                                   femtov0s::StrangeTofPosDau,
                                   femtov0s::StrangeTofNegDau);
using FXiExtras = FXiExtras_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegas_001, "FOMEGA", 1, //! femto omegas
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtocascades::BachelorId,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FOmegas = FOmegas_001;
using FOmega = FOmegas::iterator;
using StoredFOmegas = StoredFOmegas_001;

// table for basic omega information, compressed/binned kinematics
DECLARE_SOA_TABLE_STAGED_VERSIONED(FLiteOmegas_001, "FLITEOMEGA", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FLiteColId,
                                   femtobase::lite::SignedBinnedPt,
                                   femtobase::lite::BinnedEta,
                                   femtobase::lite::BinnedPhi,
                                   femtocascades::lite::BinnedOmegaMass,
                                   femtocascades::lite::BachelorId,
                                   femtov0s::lite::PosDauId,
                                   femtov0s::lite::NegDauId,
                                   femtobase::lite::Sign<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::SignedPt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::signedpt::Pt<femtobase::lite::SignedBinnedPt>,
                                   femtobase::lite::Eta<femtobase::lite::BinnedEta>,
                                   femtobase::lite::Phi<femtobase::lite::BinnedPhi>,
                                   femtocascades::lite::OmegaMass<femtocascades::lite::BinnedOmegaMass>);
using FLiteOmegas = FLiteOmegas_001;
using FLiteOmega = FLiteOmegas::iterator;
using StoredFLiteOmegas = StoredFLiteOmegas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegaMasks_001, "FOMEGAMASK", 1, //! legacy omega masks
                                   femtocascades::legacy001::Mask);

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegaMasks_002, "FOMEGAMASK", 2, //! omega masks
                                   femtocascades::Mask);
using FOmegaMasks = FOmegaMasks_002;
using StoredFOmegaMasks = StoredFOmegaMasks_002;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegaExtras_001, "FOMEGAEXTRA", 1, //! omega extra information
                                   femtocascades::MassXi,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaMass,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv,
                                   femtocascades::StrangeTofBachelor,
                                   femtov0s::StrangeTofPosDau,
                                   femtov0s::StrangeTofNegDau);
using FOmegaExtras = FOmegaExtras_001;

namespace femtocharmhadrons
{
// bitmask column
DECLARE_SOA_COLUMN(Mask, mask, o2::analysis::femto::datatypes::CharmHadronMaskType); //! selection bitmask

// daughter links: row indices into the femto track table
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FTracks, "_PosDau"); //! + prong (pion in D0)
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FTracks, "_NegDau"); //! - prong (kaon in D0)

// QA/debug columns
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float); //! d0*d0 of the two prongs
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
// ML BDT scores: [0] background, [1] prompt (D0 from c), [2] non-prompt (D0 from b decay)
DECLARE_SOA_COLUMN(MlProbD0Bkg, mlProbD0Bkg, float);                   //! D0 hypothesis: background score
DECLARE_SOA_COLUMN(MlProbD0Prompt, mlProbD0Prompt, float);             //! D0 hypothesis: prompt score
DECLARE_SOA_COLUMN(MlProbD0NonPrompt, mlProbD0NonPrompt, float);       //! D0 hypothesis: non-prompt score
DECLARE_SOA_COLUMN(MlProbD0barBkg, mlProbD0barBkg, float);             //! D0bar hypothesis: background
DECLARE_SOA_COLUMN(MlProbD0barPrompt, mlProbD0barPrompt, float);       //! D0bar hypothesis: prompt
DECLARE_SOA_COLUMN(MlProbD0barNonPrompt, mlProbD0barNonPrompt, float); //! D0bar hypothesis: non-prompt
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int8_t);                          //! PWGHF selection flag
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int8_t);                    //! PWGHF selection flag
} // namespace femtocharmhadrons

DECLARE_SOA_TABLE_STAGED_VERSIONED(FD0s_001, "FD0", 1, //! femto D0/D0bar (kinematics only)
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId,
                                   femtobase::stored::SignedPt, //! sign encodes D0(+)/D0bar(-)
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass, //! mass of the accepted hypothesis
                                   femtocharmhadrons::PosDauId,
                                   femtocharmhadrons::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FD0s = FD0s_001;
using StoredFD0s = StoredFD0s_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FD0Masks_001, "FD0MASK", 1, //! femto D0 selection bitmask
                                   femtocharmhadrons::Mask);
using FD0Masks = FD0Masks_001;
using StoredFD0Masks = StoredFD0Masks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FD0Extras_001, "FD0EXTRA", 1, //! femto D0 QA / debug
                                   femtocharmhadrons::Cpa,
                                   femtocharmhadrons::CpaXY,
                                   femtocharmhadrons::DecayLength,
                                   femtocharmhadrons::DecayLengthXY,
                                   femtocharmhadrons::ImpactParameterProduct,
                                   femtocharmhadrons::CosThetaStar,
                                   femtocharmhadrons::MlProbD0Bkg,
                                   femtocharmhadrons::MlProbD0Prompt,
                                   femtocharmhadrons::MlProbD0NonPrompt,
                                   femtocharmhadrons::MlProbD0barBkg,
                                   femtocharmhadrons::MlProbD0barPrompt,
                                   femtocharmhadrons::MlProbD0barNonPrompt,
                                   femtocharmhadrons::IsSelD0, // raw PWGHF selection flags
                                   femtocharmhadrons::IsSelD0bar);
using FD0Extras = FD0Extras_001;

// tables for monte carlo
namespace femtomccollisions
{
// DECLARE_SOA_COLUMN(Mult, mult, int);   //! Multiplicity of the event as given by the generator in |eta|<0.8
// DECLARE_SOA_COLUMN(Cent, cent, float);
//
} // namespace femtomccollisions

DECLARE_SOA_TABLE_STAGED_VERSIONED(FMcCols_001, "FMCCOL", 1, //! femto mc collisions
                                   o2::soa::Index<>,
                                   femtocollisions::PosZ, //! Multiplicity of the event as given by the generator in |eta|<0.8
                                   femtocollisions::Mult,
                                   femtocollisions::Cent);
using FMcCols = FMcCols_001;
using FMcCol = FMcCols_001::iterator;

namespace femtomcparticle
{
DECLARE_SOA_COLUMN(Origin, origin, o2::analysis::femto::datatypes::McOriginType); //! Multiplicity of the event as given by the generator in |eta|<0.8
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);                                        //! Multiplicity of the event as given by the generator in |eta|<0.8
DECLARE_SOA_INDEX_COLUMN(FMcCol, fMcCol);                                         //!
} // namespace femtomcparticle

// table for basic track information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FMcParticles_001, "FMCPARTICLE", 1, //! femto tracks
                                   o2::soa::Index<>,
                                   femtomcparticle::FMcColId,
                                   femtomcparticle::Origin,
                                   femtomcparticle::PdgCode,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FMcParticles = FMcParticles_001;
using FMcParticle = FMcParticles::iterator;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FMcMothers_001, "FMCMOTHER", 1, //! first direct mother of the monte carlo particle
                                   o2::soa::Index<>,               // no collision index needed since the mother is retrieved from the daughter mc particle
                                   femtomcparticle::Origin,
                                   femtomcparticle::PdgCode,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Phi>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);

using FMcMothers = FMcMothers_001;
using FMcMother = FMcMothers::iterator;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FMcPartMoths_001, "FMCPARTMOTH", 1, //! first partonic mother of the monte carlo particle after hadronization
                                   o2::soa::Index<>,
                                   femtomcparticle::PdgCode);
using FMcPartMoths = FMcPartMoths_001;
using FMcPartMoth = FMcPartMoths::iterator;

namespace femtolabels
{
DECLARE_SOA_INDEX_COLUMN(FMcCol, fMcCol);           //!
DECLARE_SOA_INDEX_COLUMN(FMcParticle, fMcParticle); //!

DECLARE_SOA_INDEX_COLUMN(FMcMother, fMcMother);     //!
DECLARE_SOA_INDEX_COLUMN(FMcPartMoth, fMcPartMoth); //!
} // namespace femtolabels

DECLARE_SOA_TABLE(FColLabels, "AOD", "FCOLMCLABEL", femtolabels::FMcColId);

DECLARE_SOA_TABLE(FTrackLabels, "AOD", "FTRACKLABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FLambdaLabels, "AOD", "FLAMBDALABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FK0shortLabels, "AOD", "FK0SHORTLABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FD0Labels, "AOD", "FD0LABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FSigmaLabels, "AOD", "FSIGMALABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FSigmaPlusLabels, "AOD", "FSIGMAPLUSLABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FXiLabels, "AOD", "FXILABEL", femtolabels::FMcParticleId);

DECLARE_SOA_TABLE(FOmegaLabels, "AOD", "FOMEGALABEL", femtolabels::FMcParticleId);

// for mc only processing, we also need Labels pointing from mc particles to mothers and partonic mothers
DECLARE_SOA_TABLE(FMcMotherLabels, "AOD", "FMCMOTHERLABEL",
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

} // namespace o2::aod
#endif // PWGCF_FEMTO_DATAMODEL_FEMTOTABLES_H_
