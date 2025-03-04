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

/// \file LFResonanceTables.h
/// \brief Definitions of tables of resonance decay candidates
///
/// Inspired by StrangenessTables.h, FemtoDerived.h
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// \author Nasir Mehdi Malik <nasir.mehdi.malik@cern.ch>
///

#ifndef PWGLF_DATAMODEL_LFRESONANCETABLES_H_
#define PWGLF_DATAMODEL_LFRESONANCETABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
enum {
  kECbegin = 0,
  kINEL = 1,
  kINEL10,
  kINELg0,
  kINELg010,
  kTrig,
  kTrig10,
  kTrigINELg0,
  kTrigINELg010,
  kSel8,
  kSel810,
  kSel8INELg0,
  kSel8INELg010,
  kAllCuts,
  kAllCuts10,
  kAllCutsINELg0,
  kAllCutsINELg010,
  kECend,
};
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, Collisions, "_Col"); //!
DECLARE_SOA_COLUMN(Cent, cent, float);                                        //! Centrality (Multiplicity) percentile (Default: FT0M)
DECLARE_SOA_COLUMN(Spherocity, spherocity, float);                            //! Spherocity of the event
DECLARE_SOA_COLUMN(EvtPl, evtPl, float);                                      //! Second harmonic event plane
DECLARE_SOA_COLUMN(EvtPlResAB, evtPlResAB, float);                            //! Second harmonic event plane resolution of A-B sub events
DECLARE_SOA_COLUMN(EvtPlResAC, evtPlResAC, float);                            //! Second harmonic event plane resolution of A-C sub events
DECLARE_SOA_COLUMN(EvtPlResBC, evtPlResBC, float);                            //! Second harmonic event plane resolution of B-C sub events
DECLARE_SOA_COLUMN(BMagField, bMagField, float);                              //! Magnetic field
// MC
DECLARE_SOA_COLUMN(IsVtxIn10, isVtxIn10, bool);               //! Vtx10
DECLARE_SOA_COLUMN(IsINELgt0, isINELgt0, bool);               //! INEL>0
DECLARE_SOA_COLUMN(IsTriggerTVX, isTriggerTVX, bool);         //! TriggerTVX
DECLARE_SOA_COLUMN(IsInSel8, isInSel8, bool);                 //! InSel8
DECLARE_SOA_COLUMN(IsInAfterAllCuts, isInAfterAllCuts, bool); //! InAfterAllCuts
DECLARE_SOA_COLUMN(ImpactParameter, impactParameter, float);  //! ImpactParameter

} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOLLISIONS",
                  o2::soa::Index<>,
                  resocollision::CollisionId,
                  o2::aod::mult::MultNTracksPV,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::Cent,
                  resocollision::Spherocity,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC,
                  resocollision::BMagField,
                  timestamp::Timestamp,
                  evsel::NumTracksInTimeRange);
using ResoCollision = ResoCollisions::iterator;

DECLARE_SOA_TABLE(ResoMCCollisions, "AOD", "RESOMCCOLLISIONS",
                  o2::soa::Index<>,
                  resocollision::IsVtxIn10,
                  resocollision::IsINELgt0,
                  resocollision::IsTriggerTVX,
                  resocollision::IsInSel8,
                  resocollision::IsInAfterAllCuts,
                  resocollision::ImpactParameter);
using ResoMCCollision = ResoMCCollisions::iterator;

DECLARE_SOA_TABLE(ResoSpheroCollisions, "AOD", "RESOSPHEROCOLLISIONS",
                  o2::soa::Index<>,
                  resocollision::CollisionId,
                  resocollision::Spherocity);
using ResoSpheroCollision = ResoSpheroCollisions::iterator;

DECLARE_SOA_TABLE(ResoEvtPlCollisions, "AOD", "RESOEVTPLCOLLISIONS",
                  o2::soa::Index<>,
                  resocollision::CollisionId,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC);
using ResoEvtPlCollision = ResoEvtPlCollisions::iterator;

// For DF mixing study
DECLARE_SOA_TABLE(ResoCollisionDFs, "AOD", "RESOCOLLISIONDFS",
                  o2::soa::Index<>,
                  // resocollision::CollisionId,
                  o2::aod::mult::MultNTracksPV,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::Cent,
                  resocollision::Spherocity,
                  resocollision::EvtPl,
                  resocollision::EvtPlResAB,
                  resocollision::EvtPlResAC,
                  resocollision::EvtPlResBC,
                  resocollision::BMagField,
                  timestamp::Timestamp,
                  evsel::NumTracksInTimeRange);
using ResoCollisionDF = ResoCollisionDFs::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_INDEX_COLUMN(ResoCollisionDF, resoCollisionDF);
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Trk");            //!
DECLARE_SOA_INDEX_COLUMN_FULL(V0, v0, int, V0s, "_V0");                      //!
DECLARE_SOA_INDEX_COLUMN_FULL(Cascade, cascade, int, Cascades, "_Cas");      //!
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                           //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                           //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                           //! p_z (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                         //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                             //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                           //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                                //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(CascadeIndices, cascadeIndices, int[3]);                  //! Field for the track indices to remove auto-correlations (ordered: positive, negative, bachelor)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                      //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);         //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);                     //! Number of TPC clusters found
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);            //! Is global track without DCA
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);                      //! Is global track
DECLARE_SOA_COLUMN(IsPrimaryTrack, isPrimaryTrack, bool);                    //! Is primary track
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);                  //! Is primary vertex contributor
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                                    //! Has TOF
DECLARE_SOA_COLUMN(TpcNSigmaPi10, tpcNSigmaPi10, int8_t);                    //! TPC PID x10 of the track as Pion
DECLARE_SOA_COLUMN(TpcNSigmaKa10, tpcNSigmaKa10, int8_t);                    //! TPC PID x10 of the track as Kaon
DECLARE_SOA_COLUMN(TpcNSigmaPr10, tpcNSigmaPr10, int8_t);                    //! TPC PID x10 of the track as Proton
DECLARE_SOA_COLUMN(TofNSigmaPi10, tofNSigmaPi10, int8_t);                    //! TOF PID x10 of the track as Pion
DECLARE_SOA_COLUMN(TofNSigmaKa10, tofNSigmaKa10, int8_t);                    //! TOF PID x10 of the track as Kaon
DECLARE_SOA_COLUMN(TofNSigmaPr10, tofNSigmaPr10, int8_t);                    //! TOF PID x10 of the track as Proton
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);                               //! DCA between daughters
DECLARE_SOA_COLUMN(CascDaughDCA, cascDaughDCA, float);                       //! DCA between daughters from cascade
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);                                 //! V0 Cosine of Pointing Angle
DECLARE_SOA_COLUMN(CascCosPA, cascCosPA, float);                             //! Cascade Cosine of Pointing Angle
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                                 //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);                         //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);                               //! The invariant mass of V0 candidate, assuming k0s
DECLARE_SOA_COLUMN(MXi, mXi, float);                                         //! The invariant mass of Xi candidate
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);                         //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(CascTransRadius, cascTransRadius, float);                 //! Transverse radius of the decay vertex from cascade
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);                             //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);                             //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);                             //! Z position of the decay vertex
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosPi, daughterTPCNSigmaPosPi, float);   //! TPC PID of the positive daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosKa, daughterTPCNSigmaPosKa, float);   //! TPC PID of the positive daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaPosPr, daughterTPCNSigmaPosPr, float);   //! TPC PID of the positive daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegPi, daughterTPCNSigmaNegPi, float);   //! TPC PID of the negative daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegKa, daughterTPCNSigmaNegKa, float);   //! TPC PID of the negative daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaNegPr, daughterTPCNSigmaNegPr, float);   //! TPC PID of the negative daughter as Proton
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachPi, daughterTPCNSigmaBachPi, float); //! TPC PID of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachKa, daughterTPCNSigmaBachKa, float); //! TPC PID of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTPCNSigmaBachPr, daughterTPCNSigmaBachPr, float); //! TPC PID of the bachelor daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosPi, daughterTOFNSigmaPosPi, float);   //! TOF PID of the positive daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosKa, daughterTOFNSigmaPosKa, float);   //! TOF PID of the positive daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaPosPr, daughterTOFNSigmaPosPr, float);   //! TOF PID of the positive daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegPi, daughterTOFNSigmaNegPi, float);   //! TOF PID of the negative daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegKa, daughterTOFNSigmaNegKa, float);   //! TOF PID of the negative daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaNegPr, daughterTOFNSigmaNegPr, float);   //! TOF PID of the negative daughter as Proton
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachPi, daughterTOFNSigmaBachPi, float); //! TOF PID of the bachelor daughter as Pion
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachKa, daughterTOFNSigmaBachKa, float); //! TOF PID of the bachelor daughter as Kaon
DECLARE_SOA_COLUMN(DaughterTOFNSigmaBachPr, daughterTOFNSigmaBachPr, float); //! TOF PID of the bachelor daughter as Proton
// For MC
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! Index of the corresponding MC particle
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MotherId, motherId, int);         //! Id of the mother particle
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int);       //! PDG code of the mother particle
DECLARE_SOA_COLUMN(DaughterPDG1, daughterPDG1, int); //! PDG code of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterPDG2, daughterPDG2, int); //! PDG code of the second Daughter particle
DECLARE_SOA_COLUMN(DaughterID1, daughterID1, int);   //! Id of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterID2, daughterID2, int);   //! Id of the second Daughter particle
DECLARE_SOA_COLUMN(SiblingIds, siblingIds, int[2]);  //! Index of the particles with the same mother
DECLARE_SOA_COLUMN(BachTrkID, bachTrkID, int);       //! Id of the bach track from cascade
DECLARE_SOA_COLUMN(V0ID, v0ID, int);                 //! Id of the V0 from cascade
// Dynamic columns
// TPC PID return value/10
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPi, tpcNSigmaPi,
                           [](int8_t tpcNSigmaPi10) {
                             return (float)tpcNSigmaPi10 / 10.f;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaKa, tpcNSigmaKa,
                           [](int8_t tpcNSigmaKa10) {
                             return (float)tpcNSigmaKa10 / 10.f;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TpcNSigmaPr, tpcNSigmaPr,
                           [](int8_t tpcNSigmaPr10) {
                             return (float)tpcNSigmaPr10 / 10.f;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPi, tofNSigmaPi,
                           [](int8_t tofNSigmaPi10) {
                             return (float)tofNSigmaPi10 / 10.f;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaKa, tofNSigmaKa,
                           [](int8_t tofNSigmaKa10) {
                             return (float)tofNSigmaKa10 / 10.f;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(TofNSigmaPr, tofNSigmaPr,
                           [](int8_t tofNSigmaPr10) {
                             return (float)tofNSigmaPr10 / 10.f;
                           });
} // namespace resodaughter
DECLARE_SOA_TABLE(ResoTracks, "AOD", "RESOTRACKS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::TrackId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  resodaughter::TPCNClsFound,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  resodaughter::HasTOF,
                  resodaughter::TpcNSigmaPi10,
                  resodaughter::TpcNSigmaKa10,
                  resodaughter::TpcNSigmaPr10,
                  resodaughter::TofNSigmaPi10,
                  resodaughter::TofNSigmaKa10,
                  resodaughter::TofNSigmaPr10,
                  o2::aod::track::TPCSignal,
                  o2::aod::track::PassedITSRefit,
                  o2::aod::track::PassedTPCRefit,
                  resodaughter::IsGlobalTrackWoDCA,
                  resodaughter::IsGlobalTrack,
                  resodaughter::IsPrimaryTrack,
                  resodaughter::IsPVContributor,
                  // Dynamic columns
                  resodaughter::TpcNSigmaPi<resodaughter::TpcNSigmaPi10>,
                  resodaughter::TpcNSigmaKa<resodaughter::TpcNSigmaKa10>,
                  resodaughter::TpcNSigmaPr<resodaughter::TpcNSigmaPr10>,
                  resodaughter::TofNSigmaPi<resodaughter::TofNSigmaPi10>,
                  resodaughter::TofNSigmaKa<resodaughter::TofNSigmaKa10>,
                  resodaughter::TofNSigmaPr<resodaughter::TofNSigmaPr10>);
using ResoTrack = ResoTracks::iterator;

// For DF mixing study
DECLARE_SOA_TABLE(ResoTrackDFs, "AOD", "RESOTRACKDFS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionDFId,
                  //  resodaughter::TrackId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  resodaughter::TPCNClsFound,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  resodaughter::HasTOF,
                  resodaughter::TpcNSigmaPi10,
                  resodaughter::TpcNSigmaKa10,
                  resodaughter::TpcNSigmaPr10,
                  resodaughter::TofNSigmaPi10,
                  resodaughter::TofNSigmaKa10,
                  resodaughter::TofNSigmaPr10,
                  o2::aod::track::TPCSignal,
                  o2::aod::track::PassedITSRefit,
                  o2::aod::track::PassedTPCRefit,
                  resodaughter::IsGlobalTrackWoDCA,
                  resodaughter::IsGlobalTrack,
                  resodaughter::IsPrimaryTrack,
                  resodaughter::IsPVContributor,
                  // Dynamic columns
                  resodaughter::TpcNSigmaPi<resodaughter::TpcNSigmaPi10>,
                  resodaughter::TpcNSigmaKa<resodaughter::TpcNSigmaKa10>,
                  resodaughter::TpcNSigmaPr<resodaughter::TpcNSigmaPr10>,
                  resodaughter::TofNSigmaPi<resodaughter::TofNSigmaPi10>,
                  resodaughter::TofNSigmaKa<resodaughter::TofNSigmaKa10>,
                  resodaughter::TofNSigmaPr<resodaughter::TofNSigmaPr10>);
using ResoTrackDF = ResoTrackDFs::iterator;

DECLARE_SOA_TABLE(ResoV0s, "AOD", "RESOV0S",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::V0Id,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Indices,
                  resodaughter::DaughterTPCNSigmaPosPi,
                  resodaughter::DaughterTPCNSigmaPosKa,
                  resodaughter::DaughterTPCNSigmaPosPr,
                  resodaughter::DaughterTPCNSigmaNegPi,
                  resodaughter::DaughterTPCNSigmaNegKa,
                  resodaughter::DaughterTPCNSigmaNegPr,
                  resodaughter::DaughterTOFNSigmaPosPi,
                  resodaughter::DaughterTOFNSigmaPosKa,
                  resodaughter::DaughterTOFNSigmaPosPr,
                  resodaughter::DaughterTOFNSigmaNegPi,
                  resodaughter::DaughterTOFNSigmaNegKa,
                  resodaughter::DaughterTOFNSigmaNegPr,
                  resodaughter::V0CosPA,
                  resodaughter::DaughDCA,
                  v0data::DCAPosToPV,
                  v0data::DCANegToPV,
                  v0data::DCAV0ToPV,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::MK0Short,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoV0 = ResoV0s::iterator;

DECLARE_SOA_TABLE(ResoCascades, "AOD", "RESOCASCADES",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::CascadeId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::CascadeIndices,
                  resodaughter::DaughterTPCNSigmaPosPi,
                  resodaughter::DaughterTPCNSigmaPosKa,
                  resodaughter::DaughterTPCNSigmaPosPr,
                  resodaughter::DaughterTPCNSigmaNegPi,
                  resodaughter::DaughterTPCNSigmaNegKa,
                  resodaughter::DaughterTPCNSigmaNegPr,
                  resodaughter::DaughterTPCNSigmaBachPi,
                  resodaughter::DaughterTPCNSigmaBachKa,
                  resodaughter::DaughterTPCNSigmaBachPr,
                  resodaughter::DaughterTOFNSigmaPosPi,
                  resodaughter::DaughterTOFNSigmaPosKa,
                  resodaughter::DaughterTOFNSigmaPosPr,
                  resodaughter::DaughterTOFNSigmaNegPi,
                  resodaughter::DaughterTOFNSigmaNegKa,
                  resodaughter::DaughterTOFNSigmaNegPr,
                  resodaughter::DaughterTOFNSigmaBachPi,
                  resodaughter::DaughterTOFNSigmaBachKa,
                  resodaughter::DaughterTOFNSigmaBachPr,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  cascdata::DCAPosToPV,
                  cascdata::DCANegToPV,
                  cascdata::DCABachToPV,
                  v0data::DCAV0ToPV,
                  cascdata::DCAXYCascToPV,
                  cascdata::DCAZCascToPV,
                  cascdata::Sign,
                  resodaughter::MLambda,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoCascade = ResoCascades::iterator;

DECLARE_SOA_TABLE(ResoCascadeDFs, "AOD", "RESOCASCADEDFS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionDFId,
                  // resodaughter::CascadeId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::CascadeIndices,
                  resodaughter::DaughterTPCNSigmaPosPi,
                  resodaughter::DaughterTPCNSigmaPosKa,
                  resodaughter::DaughterTPCNSigmaPosPr,
                  resodaughter::DaughterTPCNSigmaNegPi,
                  resodaughter::DaughterTPCNSigmaNegKa,
                  resodaughter::DaughterTPCNSigmaNegPr,
                  resodaughter::DaughterTPCNSigmaBachPi,
                  resodaughter::DaughterTPCNSigmaBachKa,
                  resodaughter::DaughterTPCNSigmaBachPr,
                  resodaughter::DaughterTOFNSigmaPosPi,
                  resodaughter::DaughterTOFNSigmaPosKa,
                  resodaughter::DaughterTOFNSigmaPosPr,
                  resodaughter::DaughterTOFNSigmaNegPi,
                  resodaughter::DaughterTOFNSigmaNegKa,
                  resodaughter::DaughterTOFNSigmaNegPr,
                  resodaughter::DaughterTOFNSigmaBachPi,
                  resodaughter::DaughterTOFNSigmaBachKa,
                  resodaughter::DaughterTOFNSigmaBachPr,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  cascdata::DCAPosToPV,
                  cascdata::DCANegToPV,
                  cascdata::DCABachToPV,
                  v0data::DCAV0ToPV,
                  cascdata::DCAXYCascToPV,
                  cascdata::DCAZCascToPV,
                  cascdata::Sign,
                  resodaughter::MLambda,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoCascadeDF = ResoCascadeDFs::iterator;

DECLARE_SOA_TABLE(ResoMCTracks, "AOD", "RESOMCTRACKS",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::SiblingIds,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCTrack = ResoMCTracks::iterator;

DECLARE_SOA_TABLE(ResoMCV0s, "AOD", "RESOMCV0S",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::DaughterID1,
                  resodaughter::DaughterID2,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCV0 = ResoMCV0s::iterator;

DECLARE_SOA_TABLE(ResoMCCascades, "AOD", "RESOMCCASCADES",
                  mcparticle::PdgCode,
                  resodaughter::MotherId,
                  resodaughter::MotherPDG,
                  resodaughter::BachTrkID,
                  resodaughter::V0ID,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCCascade = ResoMCCascades::iterator;

DECLARE_SOA_TABLE(ResoMCParents, "AOD", "RESOMCPARENTS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::McParticleId,
                  mcparticle::PdgCode,
                  resodaughter::DaughterPDG1,
                  resodaughter::DaughterPDG2,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  mcparticle::Y,
                  mcparticle::E,
                  mcparticle::StatusCode);
using ResoMCParent = ResoMCParents::iterator;

using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>; // without Extra
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>; // Without Extra

using ResoCollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
using ResoRun2CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
using ResoCollisionCandidatesMC = soa::Join<ResoCollisionCandidates, aod::McCollisionLabels>;
using ResoRun2CollisionCandidatesMC = soa::Join<ResoRun2CollisionCandidates, aod::McCollisionLabels>;
using ResoTrackCandidates = aod::Reso2TracksPIDExt;
using ResoTrackCandidatesMC = soa::Join<ResoTrackCandidates, aod::McTrackLabels>;
using ResoV0Candidates = aod::V0Datas;
using ResoV0CandidatesMC = soa::Join<ResoV0Candidates, aod::McV0Labels>;
using ResoCascadesCandidates = aod::CascDatas;
using ResoCascadesCandidatesMC = soa::Join<ResoCascadesCandidates, aod::McCascLabels>;
using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFRESONANCETABLES_H_
