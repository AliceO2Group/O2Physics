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
/// \file LFOmega2012MlTables.h
/// \brief Derived Omega(2012) training and audit tables, separate for the XiK0s and Xi1530K decay modes
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// Mode A (XiK0s): Omega(2012)- -> Xi- K0S. Mode B (Xi1530K): Omega(2012)- -> Xi(1530)0 K- -> Xi- pi+ K-.
/// Only the input-object tables (events, cascades, V0s, tracks) are common; every event row belongs to one
/// decay mode (column omDecayMode), so the object rows of the two modes never mix either.
///
/// Cumulative pass bits (each bit requires all lower bits):
///  mode A: 1 loose (valid canonical candidate inside the rapidity window), 2 Xi selection, 4 K0s selection,
///          8 Xi-K0s kinematic (opening-angle) cut; selected = 15.
///  mode B: 1 loose (valid canonical candidate inside the rapidity window), 2 Xi selection, 4 track quality of the
///          pion and the kaon, 8 TOF requirement and PID of the pion and the kaon, 16 Xi(1530) mass window;
///          selected = 31.
#ifndef PWGLF_DATAMODEL_LFOMEGA2012MLTABLES_H_
#define PWGLF_DATAMODEL_LFOMEGA2012MLTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace omega2012ml
{
// Persisted relations target only these derived tables; source AO2D row numbers are scalar audit values.
// Events
DECLARE_SOA_COLUMN(OmDecayMode, omDecayMode, uint8_t);             //! 0 XiK0s, 1 Xi1530K: the mode whose candidates reference this event row
DECLARE_SOA_COLUMN(OmRecoCollisionId, omRecoCollisionId, int64_t); //! row of the reduced ResoCollisions_001 collision in the input DF
DECLARE_SOA_COLUMN(OmPosZ, omPosZ, float);
DECLARE_SOA_COLUMN(OmBField, omBField, float);
DECLARE_SOA_COLUMN(OmCentrality, omCentrality, float);
DECLARE_SOA_COLUMN(OmMultiplicity, omMultiplicity, float);
DECLARE_SOA_COLUMN(OmRecINELgt0, omRecINELgt0, bool);
// Common object columns
DECLARE_SOA_COLUMN(OmSourceRow, omSourceRow, int64_t); //! row of the object in the input ResoCascades / ResoV0s table
DECLARE_SOA_COLUMN(OmPx, omPx, float);
DECLARE_SOA_COLUMN(OmPy, omPy, float);
DECLARE_SOA_COLUMN(OmPz, omPz, float);
DECLARE_SOA_COLUMN(OmV0CosPA, omV0CosPA, float);
DECLARE_SOA_COLUMN(OmV0DaughDCA, omV0DaughDCA, float);
DECLARE_SOA_COLUMN(OmDcaPosToPV, omDcaPosToPV, float);
DECLARE_SOA_COLUMN(OmDcaNegToPV, omDcaNegToPV, float);
DECLARE_SOA_COLUMN(OmDcaV0ToPV, omDcaV0ToPV, float);
DECLARE_SOA_COLUMN(OmV0Radius, omV0Radius, float);
DECLARE_SOA_COLUMN(OmMassLambda, omMassLambda, float);
DECLARE_SOA_COLUMN(OmDecayVtxX, omDecayVtxX, float);
DECLARE_SOA_COLUMN(OmDecayVtxY, omDecayVtxY, float);
DECLARE_SOA_COLUMN(OmDecayVtxZ, omDecayVtxZ, float);
DECLARE_SOA_COLUMN(OmTpcPosPi10, omTpcPosPi10, int8_t); //! TPC nSigma x10 of the positive daughter (pion hypothesis), as in the reduced tables
DECLARE_SOA_COLUMN(OmTpcPosKa10, omTpcPosKa10, int8_t);
DECLARE_SOA_COLUMN(OmTpcPosPr10, omTpcPosPr10, int8_t);
DECLARE_SOA_COLUMN(OmTpcNegPi10, omTpcNegPi10, int8_t);
DECLARE_SOA_COLUMN(OmTpcNegKa10, omTpcNegKa10, int8_t);
DECLARE_SOA_COLUMN(OmTpcNegPr10, omTpcNegPr10, int8_t);
DECLARE_SOA_COLUMN(OmTofPosPi10, omTofPosPi10, int8_t); //! TOF nSigma x10 of the positive daughter (pion hypothesis), as in the reduced tables
DECLARE_SOA_COLUMN(OmTofPosKa10, omTofPosKa10, int8_t);
DECLARE_SOA_COLUMN(OmTofPosPr10, omTofPosPr10, int8_t);
DECLARE_SOA_COLUMN(OmTofNegPi10, omTofNegPi10, int8_t);
DECLARE_SOA_COLUMN(OmTofNegKa10, omTofNegKa10, int8_t);
DECLARE_SOA_COLUMN(OmTofNegPr10, omTofNegPr10, int8_t);
DECLARE_SOA_COLUMN(OmCrossedRowsPos, omCrossedRowsPos, uint8_t);
DECLARE_SOA_COLUMN(OmCrossedRowsNeg, omCrossedRowsNeg, uint8_t);
// Cascade-only columns
DECLARE_SOA_COLUMN(OmSign, omSign, int8_t);
DECLARE_SOA_COLUMN(OmMassXi, omMassXi, float);
DECLARE_SOA_COLUMN(OmCascCosPA, omCascCosPA, float);
DECLARE_SOA_COLUMN(OmCascDaughDCA, omCascDaughDCA, float);
DECLARE_SOA_COLUMN(OmDcaBachToPV, omDcaBachToPV, float);
DECLARE_SOA_COLUMN(OmDcaXYCascToPV, omDcaXYCascToPV, float);
DECLARE_SOA_COLUMN(OmDcaZCascToPV, omDcaZCascToPV, float);
DECLARE_SOA_COLUMN(OmCascRadius, omCascRadius, float);
DECLARE_SOA_COLUMN(OmTpcBachPi10, omTpcBachPi10, int8_t);
DECLARE_SOA_COLUMN(OmTpcBachKa10, omTpcBachKa10, int8_t);
DECLARE_SOA_COLUMN(OmTpcBachPr10, omTpcBachPr10, int8_t);
DECLARE_SOA_COLUMN(OmTofBachPi10, omTofBachPi10, int8_t);
DECLARE_SOA_COLUMN(OmTofBachKa10, omTofBachKa10, int8_t);
DECLARE_SOA_COLUMN(OmTofBachPr10, omTofBachPr10, int8_t);
DECLARE_SOA_COLUMN(OmCrossedRowsBach, omCrossedRowsBach, uint8_t);
// The daughter-ID arrays mirror the persistent ResoCascades / ResoV0s columns.
DECLARE_SOA_COLUMN(OmCascadeIndices, omCascadeIndices, int[3]); //! source track IDs of the cascade daughters (positive, negative, bachelor)
// V0-only columns
DECLARE_SOA_COLUMN(OmMassK0Short, omMassK0Short, float);
DECLARE_SOA_COLUMN(OmMassAntiLambda, omMassAntiLambda, float);
DECLARE_SOA_COLUMN(OmArmAlpha, omArmAlpha, float);
DECLARE_SOA_COLUMN(OmArmQt, omArmQt, float);
DECLARE_SOA_COLUMN(OmV0Indices, omV0Indices, int[2]); //! source track IDs of the V0 daughters (positive, negative)
// Track columns
DECLARE_SOA_COLUMN(OmSourceTrackId, omSourceTrackId, int64_t); //! trackId of the reduced micro track
DECLARE_SOA_COLUMN(OmPidPi, omPidPi, uint8_t);
DECLARE_SOA_COLUMN(OmPidKa, omPidKa, uint8_t);
DECLARE_SOA_COLUMN(OmPidPr, omPidPr, uint8_t);
DECLARE_SOA_COLUMN(OmSelectionFlags, omSelectionFlags, uint8_t);
DECLARE_SOA_COLUMN(OmTrackFlags, omTrackFlags, uint8_t);
DECLARE_SOA_COLUMN(OmCrossedRows, omCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(OmItsClusterMap, omItsClusterMap, uint8_t);
// Candidate columns
DECLARE_SOA_COLUMN(OmMass, omMass, float);
DECLARE_SOA_COLUMN(OmPt, omPt, float);
DECLARE_SOA_COLUMN(OmY, omY, float);
DECLARE_SOA_COLUMN(OmEta, omEta, float);
DECLARE_SOA_COLUMN(OmPhi, omPhi, float);
DECLARE_SOA_COLUMN(OmOpeningAngle, omOpeningAngle, float); //! Xi-K0s opening angle alpha_oa of the kinematic cut
DECLARE_SOA_COLUMN(OmCharge, omCharge, int8_t);            //! sign of the Xi
DECLARE_SOA_COLUMN(OmMassXiPi, omMassXiPi, float);
DECLARE_SOA_COLUMN(OmMassXiK, omMassXiK, float);
DECLARE_SOA_COLUMN(OmMassPiK, omMassPiK, float);
DECLARE_SOA_COLUMN(OmChargePattern, omChargePattern, uint8_t);        //! 0 signal (pion opposite to the Xi, kaon equal), 1 wrong-sign pion, 2 wrong-sign kaon, 3 both
DECLARE_SOA_COLUMN(OmPassBits, omPassBits, uint16_t);                 //! cumulative pass bits of the mode, see the file header
DECLARE_SOA_COLUMN(OmXiK0sFeatures, omXiK0sFeatures, float[63]);      //! master features of the XiK0s contract (PWGLF/Core/Omega2012MlFeatures.h)
DECLARE_SOA_COLUMN(OmXi1530KFeatures, omXi1530KFeatures, float[126]); //! master features of the Xi1530K contract (PWGLF/Core/Omega2012MlFeatures.h)
DECLARE_SOA_COLUMN(OmFeatureStatus, omFeatureStatus, uint8_t);        //! o2::analysis::omega2012ml::BuildStatus; only Ok (0) rows are written
// Truth and generated audit
DECLARE_SOA_COLUMN(OmTruthStatus, omTruthStatus, uint8_t); //! 0 data, 1 matched, 2 unmatched
DECLARE_SOA_COLUMN(OmMotherPdg, omMotherPdg, int32_t);
DECLARE_SOA_COLUMN(OmMotherId, omMotherId, int64_t);
DECLARE_SOA_COLUMN(OmXiMotherPdg, omXiMotherPdg, int32_t); //! immediate-mother PDG of the Xi (any candidate)
DECLARE_SOA_COLUMN(OmV0MotherPdg, omV0MotherPdg, int32_t); //! immediate-mother PDG of the V0 (any candidate; 311 shows a K0bar intermediate)
DECLARE_SOA_COLUMN(OmXi1530Id, omXi1530Id, int64_t);       //! immediate-mother ID of the Xi (the Xi(1530)0 candidate)
DECLARE_SOA_COLUMN(OmOriginalMcParticleId, omOriginalMcParticleId, int64_t);
DECLARE_SOA_COLUMN(OmPdg, omPdg, int32_t);
DECLARE_SOA_COLUMN(OmDaughterPdg1, omDaughterPdg1, int32_t);
DECLARE_SOA_COLUMN(OmDaughterPdg2, omDaughterPdg2, int32_t);
DECLARE_SOA_COLUMN(OmGenChannel, omGenChannel, uint8_t); //! immediate channel of a generated Omega(2012): 0 other, 1 XiK0s, 2 Xi1530K
DECLARE_SOA_COLUMN(OmGenPt, omGenPt, float);
DECLARE_SOA_COLUMN(OmGenY, omGenY, float);
} // namespace omega2012ml

DECLARE_SOA_TABLE(Omega2012MlEvents, "AOD", "OMMLEVENT",
                  o2::soa::Index<>, omega2012ml::OmDecayMode, omega2012ml::OmRecoCollisionId,
                  omega2012ml::OmPosZ, omega2012ml::OmBField, omega2012ml::OmCentrality,
                  omega2012ml::OmMultiplicity, omega2012ml::OmRecINELgt0);
namespace omega2012ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlEvent, omega2012MlEvent, int, Omega2012MlEvents, "");
} // namespace omega2012ml

DECLARE_SOA_TABLE(Omega2012MlCascades, "AOD", "OMMLCASCADE",
                  o2::soa::Index<>, omega2012ml::Omega2012MlEventId, omega2012ml::OmSourceRow,
                  omega2012ml::OmPx, omega2012ml::OmPy, omega2012ml::OmPz, omega2012ml::OmSign,
                  omega2012ml::OmMassXi, omega2012ml::OmMassLambda,
                  omega2012ml::OmV0CosPA, omega2012ml::OmCascCosPA, omega2012ml::OmV0DaughDCA, omega2012ml::OmCascDaughDCA,
                  omega2012ml::OmDcaPosToPV, omega2012ml::OmDcaNegToPV, omega2012ml::OmDcaBachToPV, omega2012ml::OmDcaV0ToPV,
                  omega2012ml::OmDcaXYCascToPV, omega2012ml::OmDcaZCascToPV,
                  omega2012ml::OmV0Radius, omega2012ml::OmCascRadius,
                  omega2012ml::OmDecayVtxX, omega2012ml::OmDecayVtxY, omega2012ml::OmDecayVtxZ,
                  omega2012ml::OmTpcPosPi10, omega2012ml::OmTpcPosKa10, omega2012ml::OmTpcPosPr10,
                  omega2012ml::OmTpcNegPi10, omega2012ml::OmTpcNegKa10, omega2012ml::OmTpcNegPr10,
                  omega2012ml::OmTpcBachPi10, omega2012ml::OmTpcBachKa10, omega2012ml::OmTpcBachPr10,
                  omega2012ml::OmTofPosPi10, omega2012ml::OmTofPosKa10, omega2012ml::OmTofPosPr10,
                  omega2012ml::OmTofNegPi10, omega2012ml::OmTofNegKa10, omega2012ml::OmTofNegPr10,
                  omega2012ml::OmTofBachPi10, omega2012ml::OmTofBachKa10, omega2012ml::OmTofBachPr10,
                  omega2012ml::OmCrossedRowsPos, omega2012ml::OmCrossedRowsNeg, omega2012ml::OmCrossedRowsBach,
                  omega2012ml::OmCascadeIndices);
DECLARE_SOA_TABLE(Omega2012MlV0s, "AOD", "OMMLV0",
                  o2::soa::Index<>, omega2012ml::Omega2012MlEventId, omega2012ml::OmSourceRow,
                  omega2012ml::OmPx, omega2012ml::OmPy, omega2012ml::OmPz,
                  omega2012ml::OmMassK0Short, omega2012ml::OmMassLambda, omega2012ml::OmMassAntiLambda,
                  omega2012ml::OmV0CosPA, omega2012ml::OmV0DaughDCA,
                  omega2012ml::OmDcaPosToPV, omega2012ml::OmDcaNegToPV, omega2012ml::OmDcaV0ToPV, omega2012ml::OmV0Radius,
                  omega2012ml::OmDecayVtxX, omega2012ml::OmDecayVtxY, omega2012ml::OmDecayVtxZ,
                  omega2012ml::OmArmAlpha, omega2012ml::OmArmQt,
                  omega2012ml::OmTpcPosPi10, omega2012ml::OmTpcPosKa10, omega2012ml::OmTpcPosPr10,
                  omega2012ml::OmTpcNegPi10, omega2012ml::OmTpcNegKa10, omega2012ml::OmTpcNegPr10,
                  omega2012ml::OmTofPosPi10, omega2012ml::OmTofPosKa10, omega2012ml::OmTofPosPr10,
                  omega2012ml::OmTofNegPi10, omega2012ml::OmTofNegKa10, omega2012ml::OmTofNegPr10,
                  omega2012ml::OmCrossedRowsPos, omega2012ml::OmCrossedRowsNeg,
                  omega2012ml::OmV0Indices);
DECLARE_SOA_TABLE(Omega2012MlTracks, "AOD", "OMMLTRACK",
                  o2::soa::Index<>, omega2012ml::Omega2012MlEventId, omega2012ml::OmSourceTrackId,
                  omega2012ml::OmPx, omega2012ml::OmPy, omega2012ml::OmPz,
                  omega2012ml::OmPidPi, omega2012ml::OmPidKa, omega2012ml::OmPidPr,
                  omega2012ml::OmSelectionFlags, omega2012ml::OmTrackFlags,
                  omega2012ml::OmCrossedRows, omega2012ml::OmItsClusterMap);
namespace omega2012ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlCascade, omega2012MlCascade, int, Omega2012MlCascades, "");
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlV0, omega2012MlV0, int, Omega2012MlV0s, "");
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlPionTrack, omega2012MlPionTrack, int, Omega2012MlTracks, "_Pion");
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlKaonTrack, omega2012MlKaonTrack, int, Omega2012MlTracks, "_Kaon");
} // namespace omega2012ml

// Mode A: Xi K0S
DECLARE_SOA_TABLE(Omega2012MlXiK0sCandidates, "AOD", "OMMLXIK0SCAND",
                  o2::soa::Index<>, omega2012ml::Omega2012MlEventId,
                  omega2012ml::Omega2012MlCascadeId, omega2012ml::Omega2012MlV0Id,
                  omega2012ml::OmMass, omega2012ml::OmPt, omega2012ml::OmY, omega2012ml::OmEta, omega2012ml::OmPhi,
                  omega2012ml::OmOpeningAngle, omega2012ml::OmCharge, omega2012ml::OmPassBits);
namespace omega2012ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlXiK0sCandidate, omega2012MlXiK0sCandidate, int, Omega2012MlXiK0sCandidates, "");
} // namespace omega2012ml
DECLARE_SOA_TABLE(Omega2012MlXiK0sInputs, "AOD", "OMMLXIK0SINPUT",
                  o2::soa::Index<>, omega2012ml::Omega2012MlXiK0sCandidateId,
                  omega2012ml::OmXiK0sFeatures, omega2012ml::OmFeatureStatus);
DECLARE_SOA_TABLE(Omega2012MlXiK0sTruth, "AOD", "OMMLXIK0STRUTH",
                  o2::soa::Index<>, omega2012ml::Omega2012MlXiK0sCandidateId,
                  omega2012ml::OmTruthStatus, omega2012ml::OmMotherPdg, omega2012ml::OmMotherId,
                  omega2012ml::OmV0MotherPdg, omega2012ml::OmXiMotherPdg);

// Mode B: Xi(1530)0 K- -> Xi- pi+ K-
DECLARE_SOA_TABLE(Omega2012MlXi1530KCandidates, "AOD", "OMMLXI1530KCAND",
                  o2::soa::Index<>, omega2012ml::Omega2012MlEventId,
                  omega2012ml::Omega2012MlCascadeId, omega2012ml::Omega2012MlPionTrackId, omega2012ml::Omega2012MlKaonTrackId,
                  omega2012ml::OmMass, omega2012ml::OmMassXiPi, omega2012ml::OmMassXiK, omega2012ml::OmMassPiK,
                  omega2012ml::OmPt, omega2012ml::OmY, omega2012ml::OmEta, omega2012ml::OmPhi,
                  omega2012ml::OmChargePattern, omega2012ml::OmPassBits);
namespace omega2012ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(Omega2012MlXi1530KCandidate, omega2012MlXi1530KCandidate, int, Omega2012MlXi1530KCandidates, "");
} // namespace omega2012ml
DECLARE_SOA_TABLE(Omega2012MlXi1530KInputs, "AOD", "OMMLXI1530KINP",
                  o2::soa::Index<>, omega2012ml::Omega2012MlXi1530KCandidateId,
                  omega2012ml::OmXi1530KFeatures, omega2012ml::OmFeatureStatus);
DECLARE_SOA_TABLE(Omega2012MlXi1530KTruth, "AOD", "OMMLXI1530KTRU",
                  o2::soa::Index<>, omega2012ml::Omega2012MlXi1530KCandidateId,
                  omega2012ml::OmTruthStatus, omega2012ml::OmMotherPdg, omega2012ml::OmMotherId,
                  omega2012ml::OmXi1530Id);

// Generated Omega(2012) parents of selected reconstructed MC events (ResoMCParents_001)
DECLARE_SOA_TABLE(Omega2012MlGenAudit, "AOD", "OMMLGENAUDIT",
                  o2::soa::Index<>, omega2012ml::OmRecoCollisionId, omega2012ml::OmOriginalMcParticleId,
                  omega2012ml::OmPdg, omega2012ml::OmDaughterPdg1, omega2012ml::OmDaughterPdg2,
                  omega2012ml::OmGenChannel, omega2012ml::OmGenPt, omega2012ml::OmGenY);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFOMEGA2012MLTABLES_H_
