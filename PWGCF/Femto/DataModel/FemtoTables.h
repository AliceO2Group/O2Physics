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

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace femtocollisions
{
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::CollisionMaskType);                //! Bitmask for collision selections
DECLARE_SOA_COLUMN(CollisionTag, collisionTag, femtodatatypes::CollisionTagType); //! Bitmask for collision selections

DECLARE_SOA_COLUMN(PosX, posX, float);             //! x coordinate of vertex
DECLARE_SOA_COLUMN(PosY, posY, float);             //! y coordinate of vertex
DECLARE_SOA_COLUMN(PosZ, posZ, float);             //! z coordinate of vertex
DECLARE_SOA_COLUMN(Mult, mult, float);             //! Multiplicity estimator set by producer
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality (~= multiplicity percentile) estimator set by producer
DECLARE_SOA_COLUMN(MagField, magField, int8_t);    //! Magnetic field in kG (5 kG at normal configuration and 2kG in low B field configuration)
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
DECLARE_SOA_COLUMN(Qn, qn, float);                 //! qn bins for dividing eventsfemtab
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

// table for collisions selections
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColMasks_001, "FCOLMASK", 1, //! track masks
                                   femtocollisions::Mask);
using FColMasks = FColMasks_001;
using StoredFColMasks = StoredFColMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FColSphericities_001, "FCOLSPHERICITY", 1, //! sphericity
                                   femtocollisions::Sphericity);
using FColSphericities = FColSphericities_001;

// table for qn values
DECLARE_SOA_TABLE_STAGED_VERSIONED(FColQns_001, "FCOLQN", 1, //! qn vector
                                   femtocollisions::Qn);
using FColQns = FColQns_001;

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
                                   cent::CentFT0C);              //! centrality from FT0C
using FColCents = FColCents_001;

namespace femtobase
// all "basic" information to perform femto analysis, i.e. collision index and kinematics
// split kinematics in stored, i.e. stored in derived data, and dynmaic, i.e. can be computed on the fly
{
namespace stored
{
// static columns
DECLARE_SOA_INDEX_COLUMN(FCol, fCol);          //! collision index of femto collision table
DECLARE_SOA_COLUMN(SignedPt, signedPt, float); //! signed pt
DECLARE_SOA_COLUMN(Pt, pt, float);             //! pt
DECLARE_SOA_COLUMN(Eta, eta, float);           //! eta
DECLARE_SOA_COLUMN(Phi, phi, float);           //! phi
DECLARE_SOA_COLUMN(Mass, mass, float);         //! mass of particle
DECLARE_SOA_COLUMN(MassAnti, massAnti, float); //! mass of antiparticle
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
// use fabs for pt so it can also be used with signed pt
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
} // namespace femtobase

namespace femtotracks
{
// columns for track selections
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::TrackMaskType); //! Bitmask for track selections

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
DECLARE_SOA_COLUMN(TofMass, tofMass, float); //! Tof mass

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

// table for track selections and PID selections
DECLARE_SOA_TABLE_STAGED_VERSIONED(FTrackMasks_001, "FTRACKMASK", 1, //! track masks
                                   femtotracks::Mask);
using FTrackMasks = FTrackMasks_001;
using StoredFTrackMasks = StoredFTrackMasks_001;

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
                                   femtotracks::TofMass,
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
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::TwoTrackResonanceMaskType); //! Bitmask for resonance selections

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
                                   femtotwotrackresonances::Mask);
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
                                   femtotwotrackresonances::Mask);
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
                                   femtotwotrackresonances::Mask);
using FRho0Masks = FRho0Masks_001;

namespace femtov0s
{
// columns for bit masks
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::V0MaskType); //! Bitmask for v0 selections

// columns for debug information
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);         //! Mass of Lambda
DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float); //! Mass of AntiLambda
DECLARE_SOA_COLUMN(MassK0short, massK0short, float);       //! Mass of K0short
DECLARE_SOA_COLUMN(CosPa, cosPa, float);                   //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(DauDca, dauDca, float);                 //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);       //! Lambda transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);           //! x coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);           //! y coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);           //! z coordinate of Lambda decay vertex
DECLARE_SOA_DYNAMIC_COLUMN(DecayVtx, decayVtx,             //! distance of decay vertex from nominal interaction point
                           [](float vtxX, float vtxY, float vtxZ) -> float {
                             return std::hypot(vtxX, vtxY, vtxZ);
                           });

// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FTracks, "_PosDau"); //! index column for positive daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FTracks, "_NegDau"); //! index column for negative daughter track

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

DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdaMasks_001, "FLAMBDAMASK", 1, //! lambda masks
                                   femtov0s::Mask);
using FLambdaMasks = FLambdaMasks_001;
using StoredFLambdaMasks = StoredFLambdaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FLambdaExtras_001, "FLAMBDAEXTRA", 1, //! lambda extra information
                                   femtobase::stored::MassAnti,          // put mass of antiparticle, i.e. antilambda mass for lambdas and vice versa
                                   femtov0s::MassK0short,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
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

DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shortMasks_001, "FK0SHORTMASK", 1, //! k0short masks
                                   femtov0s::Mask);
using FK0shortMasks = FK0shortMasks_001;
using StoredFK0shortMasks = StoredFK0shortMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FK0shortExtras_001, "FK0SHORTEXTRA", 1, //! k0short extra information
                                   femtov0s::MassLambda,
                                   femtov0s::MassAntiLambda,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
                                   femtov0s::TransRadius,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::DecayVtx<femtov0s::DecayVtxX, femtov0s::DecayVtxY, femtov0s::DecayVtxZ>);

using FK0shortExtras = FK0shortExtras_001;

namespace femtokinks
{
// columns for bit masks
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::KinkMaskType); //! Bitmask for kink selections

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
} // namespace femtokinks

// table for basic sigma information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmas_001, "FSIGMA", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::FColId, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
                                   femtobase::stored::SignedPt,
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
using FSigmas = FSigmas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FSigmaMasks_001, "FSIGMAMASKS", 1,
                                   femtokinks::Mask);
using FSigmaMasks = FSigmaMasks_001;

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
                                   femtobase::stored::FColId, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
                                   femtobase::stored::SignedPt,
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
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::CascadeMaskType); //! Bitmask for cascade selections

// columns for cascad debug information
DECLARE_SOA_COLUMN(MassXi, massXi, float);                         //! Mass of xi
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);                   //! Mass of omega
DECLARE_SOA_COLUMN(CascadeCosPa, cascadeCosPa, float);             //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(CascadeDauDca, cascadeDauDca, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(CascadeTransRadius, cascadeTransRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaCosPa, lambdaCosPa, float);               //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(LambdaDauDca, lambdaDauDca, float);             //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(LambdaTransRadius, lambdaTransRadius, float);   //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaDcaToPv, lambdaDcaToPv, float);           //! Lambda transvers radius

// id columns for bachelor
// following same style as strangeness tables were we do not store the id of the lambda, but its daughters
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int32_t, FTracks, "_Bachelor"); //! bachelor id

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

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXiMasks_001, "FXIMASK", 1, //! xi masks
                                   femtocascades::Mask);
using FXiMasks = FXiMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FXiExtras_001, "FXIEXTRA", 1, //! xi extra information
                                   femtocascades::MassOmega,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv);
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

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegaMasks_001, "FOMEGAMASK", 1, //! omega masks
                                   femtocascades::Mask);
using FOmegaMasks = FOmegaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FOmegaExtras_001, "FOMEGAEXTRA", 1, //! omega extra information
                                   femtocascades::MassXi,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv);
using FOmegaExtras = FOmegaExtras_001;

// tables for monte carlo

namespace femtomccollisions
{
DECLARE_SOA_COLUMN(MultMc, multMc, int);   //! Multiplicity of the event as given by the generator in |eta|<0.8
DECLARE_SOA_COLUMN(CentMc, centMc, float); //! Multiplicity of the event as given by the generator in |eta|<0.8
                                           //
} // namespace femtomccollisions

DECLARE_SOA_TABLE_STAGED_VERSIONED(FMcCols_001, "FMCCOL", 1, //! femto mc collisions
                                   o2::soa::Index<>,
                                   femtomccollisions::MultMc,
                                   femtomccollisions::CentMc);
using FMcCols = FMcCols_001;
using FMcCol = FMcCols_001::iterator;

namespace femtomcparticle
{
DECLARE_SOA_COLUMN(Origin, origin, femtodatatypes::McOriginType); //! Multiplicity of the event as given by the generator in |eta|<0.8
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);                        //! Multiplicity of the event as given by the generator in |eta|<0.8
DECLARE_SOA_INDEX_COLUMN(FMcCol, fMcCol);                         //!
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
                                   o2::soa::Index<>,
                                   femtomcparticle::PdgCode);
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

DECLARE_SOA_TABLE(FColLabels, "AOD", "FCOLMCLABEL",
                  femtolabels::FMcColId);

DECLARE_SOA_TABLE(FTrackLabels, "AOD", "FTRACKLABEL",
                  femtolabels::FMcParticleId,
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

DECLARE_SOA_TABLE(FLambdaLabels, "AOD", "FLAMBDALABEL",
                  femtolabels::FMcParticleId,
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

DECLARE_SOA_TABLE(FK0shortLabels, "AOD", "FK0SHORTLABEL",
                  femtolabels::FMcParticleId,
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

DECLARE_SOA_TABLE(FSigmaLabels, "AOD", "FSIGMALABEL",
                  femtolabels::FMcParticleId,
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

DECLARE_SOA_TABLE(FSigmaPlusLabels, "AOD", "FSIGMAPLUSLABEL",
                  femtolabels::FMcParticleId,
                  femtolabels::FMcMotherId,
                  femtolabels::FMcPartMothId);

} // namespace o2::aod
#endif // PWGCF_FEMTO_DATAMODEL_FEMTOTABLES_H_
