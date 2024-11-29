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
///  Nasir Mehdi Malik

#ifndef PWGLF_DATAMODEL_LFRESONANCETABLESMERGEDF_H_
#define PWGLF_DATAMODEL_LFRESONANCETABLESMERGEDF_H_

#include <Framework/ASoA.h>
#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollisiondf
{
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality (Multiplicity) percentile (Default: FT0M)
DECLARE_SOA_COLUMN(Spherocity, spherocity, float); //! Spherocity of the event
DECLARE_SOA_COLUMN(EvtPl, evtPl, float);           //! Second harmonic event plane
DECLARE_SOA_COLUMN(EvtPlResAB, evtPlResAB, float); //! Second harmonic event plane resolution of A-B sub events
DECLARE_SOA_COLUMN(EvtPlResAC, evtPlResAC, float); //! Second harmonic event plane resolution of A-C sub events
DECLARE_SOA_COLUMN(EvtPlResBC, evtPlResBC, float); //! Second harmonic event plane resolution of B-C sub events
DECLARE_SOA_COLUMN(BMagField, bMagField, float);   //! Magnetic field
} // namespace resocollisiondf
DECLARE_SOA_TABLE(ResoCollisionDFs, "AOD", "RESOCOLLISIONDF",
                  o2::soa::Index<>,
                  o2::aod::mult::MultNTracksPV,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  resocollisiondf::Cent,
                  resocollisiondf::Spherocity,
                  resocollisiondf::EvtPl,
                  resocollisiondf::EvtPlResAB,
                  resocollisiondf::EvtPlResAC,
                  resocollisiondf::EvtPlResBC,
                  resocollisiondf::BMagField,
                  timestamp::Timestamp,
                  evsel::NumTracksInTimeRange);
using ResoCollisionDF = ResoCollisionDFs::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughterdf
{

DECLARE_SOA_INDEX_COLUMN(ResoCollisionDF, resoCollisiondf);
DECLARE_SOA_COLUMN(Pt, pt, float);                                   //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Px, px, float);                                   //! p_x (GeV/c)
DECLARE_SOA_COLUMN(Py, py, float);                                   //! p_y (GeV/c)
DECLARE_SOA_COLUMN(Pz, pz, float);                                   //! p_z (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                 //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                 //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                     //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                   //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                        //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(CascadeIndices, cascIndices, int[3]);             //! Field for the track indices to remove auto-correlations (ordered: positive, negative, bachelor)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                              //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(TPCNClsFound, tpcNClsFound, uint8_t);             //! Number of TPC clusters found
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                       //! Number of ITS clusters found
DECLARE_SOA_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, bool);    //! Is global track without DCA
DECLARE_SOA_COLUMN(IsGlobalTrack, isGlobalTrack, bool);              //! Is global track
DECLARE_SOA_COLUMN(IsPrimaryTrack, isPrimaryTrack, bool);            //! Is primary track
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);          //! Is primary vertex contributor
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool); //! Has TOF
DECLARE_SOA_COLUMN(TPCCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);               //! DCA between daughters
DECLARE_SOA_COLUMN(CascDaughDCA, cascdaughDCA, float);       //! DCA between daughters from cascade
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);                 //! V0 Cosine of Pointing Angle
DECLARE_SOA_COLUMN(CascCosPA, cascCosPA, float);             //! Cascade Cosine of Pointing Angle
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                 //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);         //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);               //! The invariant mass of V0 candidate, assuming k0s
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
} // namespace resodaughterdf
DECLARE_SOA_TABLE(ResoTrackDFs, "AOD", "RESOTRACKDFs",
                  o2::soa::Index<>,
                  resodaughterdf::ResoCollisionDFId,
                  resodaughterdf::Pt,
                  resodaughterdf::Px,
                  resodaughterdf::Py,
                  resodaughterdf::Pz,
                  resodaughterdf::Eta,
                  resodaughterdf::Phi,
                  resodaughterdf::Sign,
                  resodaughterdf::TPCNClsCrossedRows,
                  resodaughterdf::TPCNClsFound,
                  resodaughterdf::ITSNCls,
                  o2::aod::track::DcaXY,
                  o2::aod::track::DcaZ,
                  o2::aod::track::X,
                  o2::aod::track::Alpha,
                  resodaughterdf::HasITS,
                  resodaughterdf::HasTPC,
                  resodaughterdf::HasTOF,
                  o2::aod::pidtpc::TPCNSigmaPi,
                  o2::aod::pidtpc::TPCNSigmaKa,
                  o2::aod::pidtpc::TPCNSigmaPr,
                  o2::aod::pidtpc::TPCNSigmaEl,
                  o2::aod::pidtof::TOFNSigmaPi,
                  o2::aod::pidtof::TOFNSigmaKa,
                  o2::aod::pidtof::TOFNSigmaPr,
                  o2::aod::pidtof::TOFNSigmaEl,
                  o2::aod::track::TPCSignal,
                  o2::aod::track::PassedITSRefit,
                  o2::aod::track::PassedTPCRefit,
                  resodaughterdf::IsGlobalTrackWoDCA,
                  resodaughterdf::IsGlobalTrack,
                  resodaughterdf::IsPrimaryTrack,
                  resodaughterdf::IsPVContributor,
                  resodaughterdf::TPCCrossedRowsOverFindableCls,
                  o2::aod::track::ITSChi2NCl,
                  o2::aod::track::TPCChi2NCl);
using ResoTrackDF = ResoTrackDFs::iterator;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFRESONANCETABLESMERGEDF_H_
