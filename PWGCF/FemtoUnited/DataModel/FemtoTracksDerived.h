// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoTracksDerived.h
/// \brief track tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_

#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

#include <cmath>

namespace o2::aod
{
namespace femtotracks
{
// columns for track selections
DECLARE_SOA_COLUMN(TrackMask, trackMask, femtodatatypes::TrackMaskType); //! Bitmask for track selections

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
DECLARE_SOA_COLUMN(TpcInnerParam, tpcInnerParam, bool);                      //! Momentum at inner wall of Tpc
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

// PID information
// ITS PID information
DECLARE_SOA_COLUMN(ItsNSigmaEl, itsNSigmaEl, float); //! Nsigma separation with the Its detector for electron
DECLARE_SOA_COLUMN(ItsNSigmaPi, itsNSigmaPi, float); //! Nsigma separation with the Its detector for pion
DECLARE_SOA_COLUMN(ItsNSigmaKa, itsNSigmaKa, float); //! Nsigma separation with the Its detector for kaon
DECLARE_SOA_COLUMN(ItsNSigmaPr, itsNSigmaPr, float); //! Nsigma separation with the Its detector for proton
DECLARE_SOA_COLUMN(ItsNSigmaDe, itsNSigmaDe, float); //! Nsigma separation with the Its detector for deuteron
DECLARE_SOA_COLUMN(ItsNSigmaTr, itsNSigmaTr, float); //! Nsigma separation with the Its detector for triton
DECLARE_SOA_COLUMN(ItsNSigmaHe, itsNSigmaHe, float); //! Nsigma separation with the Its detector for helium3

// TPC PID information
DECLARE_SOA_COLUMN(TpcNSigmaEl, tpcNSigmaEl, float); //! Nsigma separation with the Tpc detector for electron
DECLARE_SOA_COLUMN(TpcNSigmaPi, tpcNSigmaPi, float); //! Nsigma separation with the Tpc detector for pion
DECLARE_SOA_COLUMN(TpcNSigmaKa, tpcNSigmaKa, float); //! Nsigma separation with the Tpc detector for kaon
DECLARE_SOA_COLUMN(TpcNSigmaPr, tpcNSigmaPr, float); //! Nsigma separation with the Tpc detector for proton
DECLARE_SOA_COLUMN(TpcNSigmaDe, tpcNSigmaDe, float); //! Nsigma separation with the Tpc detector for deuteron
DECLARE_SOA_COLUMN(TpcNSigmaTr, tpcNSigmaTr, float); //! Nsigma separation with the Tpc detector for triton
DECLARE_SOA_COLUMN(TpcNSigmaHe, tpcNSigmaHe, float); //! Nsigma separation with the Tpc detector for helium3

// TOF PID information
DECLARE_SOA_COLUMN(TofNSigmaEl, tofNSigmaEl, float); //! Nsigma separation with the Tof detector for electron
DECLARE_SOA_COLUMN(TofNSigmaPi, tofNSigmaPi, float); //! Nsigma separation with the Tof detector for pion
DECLARE_SOA_COLUMN(TofNSigmaKa, tofNSigmaKa, float); //! Nsigma separation with the Tof detector for kaon
DECLARE_SOA_COLUMN(TofNSigmaPr, tofNSigmaPr, float); //! Nsigma separation with the Tof detector for proton
DECLARE_SOA_COLUMN(TofNSigmaDe, tofNSigmaDe, float); //! Nsigma separation with the Tof detector for deuteron
DECLARE_SOA_COLUMN(TofNSigmaTr, tofNSigmaTr, float); //! Nsigma separation with the Tof detector for triton
DECLARE_SOA_COLUMN(TofNSigmaHe, tofNSigmaHe, float); //! Nsigma separation with the Tof detector for helium3

DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaEl, tpctofNSigmaEl, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaPi, tpctofNSigmaPi, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaKa, tpctofNSigmaKa, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaPr, tpctofNSigmaPr, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaDe, tpctofNSigmaDe, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaTr, tpctofNSigmaTr, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!
DECLARE_SOA_DYNAMIC_COLUMN(TpctofNSigmaHe, tpctofNSigmaHe, [](float tpc, float tof) -> float { return std::hypot(tpc, tof); }); //!

} // namespace femtotracks

// table for basic track information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUTracks_001, "FUTRACKS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUTracks = FUTracks_001;

// table for track selections and PID selections
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUTrackMasks_001, "FUTRACKMASKS", 1,
                                   femtotracks::TrackMask);
using FUTrackMasks = FUTrackMasks_001;

// table for track DCA
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUTrackDCAs_001, "FUTRACKDCAS", 1,
                                   femtotracks::DcaXY,
                                   femtotracks::DcaZ,
                                   femtotracks::Dca<femtotracks::DcaXY, femtotracks::DcaZ>);
using FUTrackDCAs = FUTrackDCAs_001;

// table for extra track information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUTrackExtras_001, "FUTRACKEXTRAS", 1,
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
using FUTrackExtras = FUTrackExtras_001;

// table for extra PID information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUTrackPids_001, "FUTRACKPIDS", 1,
                                   femtotracks::ItsNSigmaEl,
                                   femtotracks::ItsNSigmaPi,
                                   femtotracks::ItsNSigmaKa,
                                   femtotracks::ItsNSigmaPr,
                                   femtotracks::ItsNSigmaDe,
                                   femtotracks::ItsNSigmaTr,
                                   femtotracks::ItsNSigmaHe,
                                   femtotracks::TpcNSigmaEl,
                                   femtotracks::TpcNSigmaPi,
                                   femtotracks::TpcNSigmaKa,
                                   femtotracks::TpcNSigmaPr,
                                   femtotracks::TpcNSigmaDe,
                                   femtotracks::TpcNSigmaTr,
                                   femtotracks::TpcNSigmaHe,
                                   femtotracks::TofNSigmaEl,
                                   femtotracks::TofNSigmaPi,
                                   femtotracks::TofNSigmaKa,
                                   femtotracks::TofNSigmaPr,
                                   femtotracks::TofNSigmaDe,
                                   femtotracks::TofNSigmaTr,
                                   femtotracks::TofNSigmaHe,
                                   femtotracks::TpctofNSigmaEl<femtotracks::TpcNSigmaEl, femtotracks::TofNSigmaEl>,
                                   femtotracks::TpctofNSigmaPi<femtotracks::TpcNSigmaPi, femtotracks::TofNSigmaPi>,
                                   femtotracks::TpctofNSigmaKa<femtotracks::TpcNSigmaKa, femtotracks::TofNSigmaKa>,
                                   femtotracks::TpctofNSigmaPr<femtotracks::TpcNSigmaPr, femtotracks::TofNSigmaPr>,
                                   femtotracks::TpctofNSigmaDe<femtotracks::TpcNSigmaDe, femtotracks::TofNSigmaDe>,
                                   femtotracks::TpctofNSigmaTr<femtotracks::TpcNSigmaTr, femtotracks::TofNSigmaTr>,
                                   femtotracks::TpctofNSigmaHe<femtotracks::TpcNSigmaHe, femtotracks::TofNSigmaHe>);

using FUTrackPids = FUTrackPids_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_
