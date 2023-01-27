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
/// \file   qaEventTrack.h
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Mattia Faggin <mattia.faggin@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief  Header file for QA tasks for the track and the event properties.
///

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

// Output table declaration
namespace o2::aod
{

namespace dpgcollision
{
DECLARE_SOA_INDEX_COLUMN(BC, bc);
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
} // namespace dpgcollision

DECLARE_SOA_TABLE(DPGCollisions, "AOD", "DPGCollisions", //! Table of the DPG collisions
                  collision::PosZ,
                  dpgcollision::IsEventReject,
                  dpgcollision::RunNumber,
                  dpgcollision::NumContrib);

namespace dpgtrack
{
DECLARE_SOA_INDEX_COLUMN(DPGCollision, dpgCollision);                                    //! Index to move from track to collision
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //! Pt
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //! Phi
DECLARE_SOA_COLUMN(PtReso, ptReso, float);                                               //! Pt resolution
DECLARE_SOA_COLUMN(Sign, sign, short);                                                   //! Sign
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);                                                //! Track has the ITS
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                                                //! Track has the TPC
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                                                //! Track has the TRD
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                                                //! Track has the TOF
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, int16_t);                                 //! Clusters found in TPC
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, int16_t);                     //! Crossed rows found in TPC
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); //! Crossed rows over findable clusters in TPC
DECLARE_SOA_COLUMN(TPCFoundOverFindableCls, tpcFoundOverFindableCls, float);             //! Found over findable clusters in TPC
DECLARE_SOA_COLUMN(TPCFractionSharedCls, tpcFractionSharedCls, float);                   //! Fraction of shared clusters in TPC
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                                           //! Clusters found in ITS
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, uint8_t);                     //! Clusters found in the inner barrel of the ITS
DECLARE_SOA_COLUMN(TOFSignal, tofSignal, float);                                         //! TOF Signal (Ev. Time subtracted)

} // namespace dpgtrack

DECLARE_SOA_TABLE(DPGTracks, "AOD", "DPGTracks", //! Table of the DPG tracks
                  dpgtrack::DPGCollisionId,
                  dpgtrack::Pt, track::TPCInnerParam, dpgtrack::Eta, dpgtrack::Phi, dpgtrack::PtReso,
                  track::Flags, dpgtrack::Sign,
                  track::DcaXY, track::DcaZ, track::Length,
                  track::ITSClusterMap,
                  track::ITSChi2NCl, track::TPCChi2NCl, track::TRDChi2, track::TOFChi2,
                  dpgtrack::HasITS, dpgtrack::HasTPC, dpgtrack::HasTRD, dpgtrack::HasTOF,
                  dpgtrack::TPCNClsFound, dpgtrack::TPCNClsCrossedRows,
                  dpgtrack::TPCCrossedRowsOverFindableCls, dpgtrack::TPCFoundOverFindableCls, dpgtrack::TPCFractionSharedCls,
                  dpgtrack::ITSNCls, dpgtrack::ITSNClsInnerBarrel, track::TPCSignal, dpgtrack::TOFSignal);

namespace dpgparticles
{
DECLARE_SOA_COLUMN(PtMC, ptMC, float);                   //! Pt MC
DECLARE_SOA_COLUMN(EtaMC, etaMC, float);                 //! Eta MC
DECLARE_SOA_COLUMN(PhiMC, phiMC, float);                 //! Phi MC
DECLARE_SOA_COLUMN(ProductionMode, productionMode, int); //! ProductionMode i.e. non matched (-1), physical primary (0), weak decay product (1) or material (2)
DECLARE_SOA_DYNAMIC_COLUMN(IsNonMatched, isNonMatched,   //! True if particle is considered a non matched particle
                           [](int mode) -> bool { return mode == -1; });
DECLARE_SOA_DYNAMIC_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, //! True if particle is considered a physical primary according to the ALICE definition
                           [](int mode) -> bool { return mode == 0; });
DECLARE_SOA_DYNAMIC_COLUMN(IsFromWeakDecay, isFromWeakDecay, //! True if particle is considered from a weak decay
                           [](int mode) -> bool { return mode == 1; });
DECLARE_SOA_DYNAMIC_COLUMN(IsFromMaterial, isFromMaterial, //! True if particle is considered from a weak decay
                           [](int mode) -> bool { return mode == 2; });

} // namespace dpgparticles

DECLARE_SOA_TABLE(DPGRecoParticles, "AOD", "DPGRecoPart", //! Table of the DPG reconstructed particles
                  dpgparticles::PtMC, dpgparticles::EtaMC, dpgparticles::PhiMC,
                  mcparticle::PdgCode, dpgparticles::ProductionMode,
                  dpgparticles::IsNonMatched<dpgparticles::ProductionMode>,
                  dpgparticles::IsPhysicalPrimary<dpgparticles::ProductionMode>,
                  dpgparticles::IsFromWeakDecay<dpgparticles::ProductionMode>,
                  dpgparticles::IsFromMaterial<dpgparticles::ProductionMode>);

DECLARE_SOA_TABLE(DPGNonRecoParticles, "AOD", "DPGNonRecoPart", //! Table of the DPG non reconstructed particles
                  dpgtrack::DPGCollisionId,
                  dpgparticles::PtMC, dpgparticles::EtaMC, dpgparticles::PhiMC,
                  mcparticle::PdgCode, dpgparticles::ProductionMode,
                  dpgparticles::IsNonMatched<dpgparticles::ProductionMode>,
                  dpgparticles::IsPhysicalPrimary<dpgparticles::ProductionMode>,
                  dpgparticles::IsFromWeakDecay<dpgparticles::ProductionMode>,
                  dpgparticles::IsFromMaterial<dpgparticles::ProductionMode>,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz);
} // namespace o2::aod
