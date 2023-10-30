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

#ifndef PWGLF_DATAMODEL_LFRESONANCETABLES_H_
#define PWGLF_DATAMODEL_LFRESONANCETABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);         //! V0M multiplicity percentile (run2: V0M, run3: FT0A/C/M)
DECLARE_SOA_COLUMN(MultFT0, multFT0, int);           //! FT0 multiplicity
DECLARE_SOA_COLUMN(Spherocity, spherocity, float);   //! Spherocity of the event
DECLARE_SOA_COLUMN(BMagField, bMagField, float);     //! Magnetic field
} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOL",
                  o2::soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollision::MultV0M,
                  resocollision::MultFT0,
                  resocollision::Spherocity,
                  resocollision::BMagField,
                  timestamp::Timestamp);
using ResoCollision = ResoCollisions::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                                     //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                     //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                     //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                     //! p_z (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                       //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                     //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                          //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);   //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);      //! Is global track without DCA
DECLARE_SOA_COLUMN(IsPrimaryTrack, isPrimaryTrack, bool);              //! Is primary track
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);            //! Is primary vertex contributor
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);                              //! Has TOF
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);               //! DCA between daughters
DECLARE_SOA_COLUMN(CascDaughDCA, cascdaughDCA, float);       //! DCA between daughters from cascade
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);                 //! V0 Cosine of Pointing Angle
DECLARE_SOA_COLUMN(CascCosPA, cascCosPA, float);             //! Cascade Cosine of Pointing Angle
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                 //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);         //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(MXi, mXi, float);                         //! The invariant mass of Xi candidate
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);         //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(CascTransRadius, casctransRadius, float); //! Transverse radius of the decay vertex from cascade
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);             //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);             //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);             //! Z position of the decay vertex
// For MC
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! Index of the corresponding MC particle
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_COLUMN(ProducedByGenerator, producedByGenerator, bool);
DECLARE_SOA_COLUMN(MothersId, motherId, int);        //! Id of the mother particle
DECLARE_SOA_COLUMN(MotherPDG, motherPDG, int);       //! PDG code of the mother particle
DECLARE_SOA_COLUMN(DaughterPDG1, daughterPDG1, int); //! PDG code of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterPDG2, daughterPDG2, int); //! PDG code of the second Daughter particle
DECLARE_SOA_COLUMN(DaughterID1, daughterId1, int);   //! Id of the first Daughter particle
DECLARE_SOA_COLUMN(DaughterID2, daughterId2, int);   //! Id of the second Daughter particle
DECLARE_SOA_COLUMN(SiblingIds, siblingIds, int[2]);  //! Index of the particles with the same mother
DECLARE_SOA_COLUMN(BachTrkID, bachtrkID, int);       //! Id of the bach track from cascade
DECLARE_SOA_COLUMN(V0ID, v0ID, int);                 //! Id of the V0 from cascade
} // namespace resodaughter
DECLARE_SOA_TABLE(ResoTracks, "AOD", "RESOTRACKS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  o2::aod::track::X,
                  o2::aod::track::Alpha,
                  resodaughter::HasTOF,
                  o2::aod::pidtpc::TPCNSigmaPi,
                  o2::aod::pidtpc::TPCNSigmaKa,
                  o2::aod::pidtpc::TPCNSigmaPr,
                  o2::aod::pidtof::TOFNSigmaPi,
                  o2::aod::pidtof::TOFNSigmaKa,
                  o2::aod::pidtof::TOFNSigmaPr,
                  o2::aod::track::TPCSignal,
                  o2::aod::track::PassedITSRefit,
                  o2::aod::track::PassedTPCRefit,
                  resodaughter::IsGlobalTrackWoDCA,
                  resodaughter::IsPrimaryTrack,
                  resodaughter::IsPVContributor,
                  resodaughter::TPCCrossedRowsOverFindableCls,
                  o2::aod::track::ITSChi2NCl,
                  o2::aod::track::TPCChi2NCl);
using ResoTrack = ResoTracks::iterator;

DECLARE_SOA_TABLE(ResoV0s, "AOD", "RESOV0S",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Indices,
                  resodaughter::V0CosPA,
                  resodaughter::DaughDCA,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoV0 = ResoV0s::iterator;

DECLARE_SOA_TABLE(ResoCascades, "AOD", "RESOCASCADES",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Px,
                  resodaughter::Py,
                  resodaughter::Pz,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::Indices,
                  resodaughter::V0CosPA,
                  resodaughter::CascCosPA,
                  resodaughter::DaughDCA,
                  resodaughter::CascDaughDCA,
                  resodaughter::MXi,
                  resodaughter::TransRadius,
                  resodaughter::CascTransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ);
using ResoCascade = ResoCascades::iterator;

DECLARE_SOA_TABLE(ResoMCTracks, "AOD", "RESOMCTRACKS",
                  mcparticle::PdgCode,
                  resodaughter::MothersId,
                  resodaughter::MotherPDG,
                  resodaughter::SiblingIds,
                  resodaughter::IsPhysicalPrimary,
                  resodaughter::ProducedByGenerator);
using ResoMCTrack = ResoMCTracks::iterator;

DECLARE_SOA_TABLE(ResoMCV0s, "AOD", "RESOMCV0S",
                  mcparticle::PdgCode,
                  resodaughter::MothersId,
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
                  resodaughter::MothersId,
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
                  mcparticle::Y);
using ResoMCParent = ResoMCParents::iterator;

using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>; // without Extra
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>; // Without Extra

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFRESONANCETABLES_H_
