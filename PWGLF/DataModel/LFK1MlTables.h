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
/// \file LFK1MlTables.h
/// \brief Derived K1 microtrack training and audit tables
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
#ifndef PWGLF_DATAMODEL_LFK1MLTABLES_H_
#define PWGLF_DATAMODEL_LFK1MLTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace k1ml
{
// Persisted relations target only these derived tables, so AO2D merge tools
// can relocate them. Source AO2D row numbers remain scalar audit values.
DECLARE_SOA_COLUMN(K1RecoCollisionId, k1RecoCollisionId, int64_t); //! row index of the reduced ResoCollisions_001 collision in the input DF
DECLARE_SOA_COLUMN(K1PosZ, k1PosZ, float);
DECLARE_SOA_COLUMN(K1BField, k1BField, float);
DECLARE_SOA_COLUMN(K1Centrality, k1Centrality, float);
DECLARE_SOA_COLUMN(K1Multiplicity, k1Multiplicity, float);
DECLARE_SOA_COLUMN(K1RecINELgt0, k1RecINELgt0, bool);
DECLARE_SOA_COLUMN(K1SourceTrackId, k1SourceTrackId, int64_t); //! trackId of the reduced micro track
DECLARE_SOA_COLUMN(K1Px, k1Px, float);
DECLARE_SOA_COLUMN(K1Py, k1Py, float);
DECLARE_SOA_COLUMN(K1Pz, k1Pz, float);
DECLARE_SOA_COLUMN(K1PidPi, k1PidPi, uint8_t);
DECLARE_SOA_COLUMN(K1PidKa, k1PidKa, uint8_t);
DECLARE_SOA_COLUMN(K1PidPr, k1PidPr, uint8_t);
DECLARE_SOA_COLUMN(K1SelectionFlags, k1SelectionFlags, uint8_t);
DECLARE_SOA_COLUMN(K1TrackFlags, k1TrackFlags, uint8_t);
DECLARE_SOA_COLUMN(K1CrossedRows, k1CrossedRows, uint8_t);
DECLARE_SOA_COLUMN(K1ItsClusterMap, k1ItsClusterMap, uint8_t);
DECLARE_SOA_COLUMN(K1Mass, k1Mass, float);
DECLARE_SOA_COLUMN(K1MassPiPi, k1MassPiPi, float);
DECLARE_SOA_COLUMN(K1MassKaPiSame, k1MassKaPiSame, float);
DECLARE_SOA_COLUMN(K1MassKaPiOpp, k1MassKaPiOpp, float);
DECLARE_SOA_COLUMN(K1ScalarSumPt, k1ScalarSumPt, float);
DECLARE_SOA_COLUMN(K1PiPiPt, k1PiPiPt, float);
DECLARE_SOA_COLUMN(K1Pt, k1Pt, float);
DECLARE_SOA_COLUMN(K1Y, k1Y, float);
DECLARE_SOA_COLUMN(K1Eta, k1Eta, float);
DECLARE_SOA_COLUMN(K1Phi, k1Phi, float);
DECLARE_SOA_COLUMN(K1Charge, k1Charge, int8_t);
// Cumulative selection bits of the unlike-sign candidate:
// 1 = valid canonical candidate inside the K1 rapidity window (loose stage),
// 2 = track quality of all three tracks,
// 4 = TOF requirement and PID of all three tracks,
// 8 = pion-pair pT and secondary mass window,
// 16 = candidate cuts.
// Candidates written at the "selected" export stage carry 31.
DECLARE_SOA_COLUMN(K1BaselinePassBits, k1BaselinePassBits, uint16_t);
DECLARE_SOA_COLUMN(K1MasterFeatures, k1MasterFeatures, float[125]); //! 125 master features of the K1 ML feature contract (PWGLF/Core/K1MlFeatures.h)
DECLARE_SOA_COLUMN(K1FeatureStatus, k1FeatureStatus, uint8_t);      //! o2::analysis::k1ml::BuildStatus; only Ok (0) rows are written
DECLARE_SOA_COLUMN(K1TruthStatus, k1TruthStatus, uint8_t);          //! 0 data, 1 matched, 2 unmatched
DECLARE_SOA_COLUMN(K1TruthChannel, k1TruthChannel, uint8_t);        //! 0 none, 1 rho K, 2 K* pi
DECLARE_SOA_COLUMN(K1MotherPdg, k1MotherPdg, int32_t);
DECLARE_SOA_COLUMN(K1MotherId, k1MotherId, int64_t);
DECLARE_SOA_COLUMN(K1GeneratedPt, k1GeneratedPt, float); //! NaN in K1MlTruth: not filled for reconstructed candidates
DECLARE_SOA_COLUMN(K1GeneratedY, k1GeneratedY, float);   //! NaN in K1MlTruth: not filled for reconstructed candidates
DECLARE_SOA_COLUMN(K1OriginalMcParticleId, k1OriginalMcParticleId, int64_t);
DECLARE_SOA_COLUMN(K1DaughterPdg1, k1DaughterPdg1, int32_t);
DECLARE_SOA_COLUMN(K1DaughterPdg2, k1DaughterPdg2, int32_t);
DECLARE_SOA_COLUMN(K1GenSelectedRecoEvent, k1GenSelectedRecoEvent, bool); //! true: parents are taken from selected reconstructed events
} // namespace k1ml

DECLARE_SOA_TABLE(K1MlEvents, "AOD", "K1MLEVENT",
                  o2::soa::Index<>, k1ml::K1RecoCollisionId, k1ml::K1PosZ, k1ml::K1BField,
                  k1ml::K1Centrality, k1ml::K1Multiplicity, k1ml::K1RecINELgt0);
namespace k1ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(K1MlEvent, k1MlEvent, int, K1MlEvents, "");
} // namespace k1ml
DECLARE_SOA_TABLE(K1MlTracks, "AOD", "K1MLTRACK",
                  o2::soa::Index<>, k1ml::K1MlEventId, k1ml::K1SourceTrackId,
                  k1ml::K1Px, k1ml::K1Py, k1ml::K1Pz,
                  k1ml::K1PidPi, k1ml::K1PidKa, k1ml::K1PidPr,
                  k1ml::K1SelectionFlags, k1ml::K1TrackFlags,
                  k1ml::K1CrossedRows, k1ml::K1ItsClusterMap);
namespace k1ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(K1MlKaonTrack, k1MlKaonTrack, int, K1MlTracks, "_Kaon");
DECLARE_SOA_INDEX_COLUMN_FULL(K1MlSamePionTrack, k1MlSamePionTrack, int, K1MlTracks, "_Same");
DECLARE_SOA_INDEX_COLUMN_FULL(K1MlOppPionTrack, k1MlOppPionTrack, int, K1MlTracks, "_Opp");
} // namespace k1ml
DECLARE_SOA_TABLE(K1MlCandidates, "AOD", "K1MLCANDIDATE",
                  o2::soa::Index<>, k1ml::K1MlEventId,
                  k1ml::K1MlKaonTrackId, k1ml::K1MlSamePionTrackId, k1ml::K1MlOppPionTrackId,
                  k1ml::K1Mass, k1ml::K1MassPiPi, k1ml::K1MassKaPiSame, k1ml::K1MassKaPiOpp,
                  k1ml::K1ScalarSumPt, k1ml::K1PiPiPt,
                  k1ml::K1Pt, k1ml::K1Y, k1ml::K1Eta, k1ml::K1Phi, k1ml::K1Charge,
                  k1ml::K1BaselinePassBits);
namespace k1ml
{
DECLARE_SOA_INDEX_COLUMN_FULL(K1MlCandidate, k1MlCandidate, int, K1MlCandidates, "");
} // namespace k1ml
DECLARE_SOA_TABLE(K1MlInputs, "AOD", "K1MLINPUT",
                  o2::soa::Index<>, k1ml::K1MlCandidateId,
                  k1ml::K1MasterFeatures, k1ml::K1FeatureStatus);
DECLARE_SOA_TABLE(K1MlTruth, "AOD", "K1MLTRUTH",
                  o2::soa::Index<>, k1ml::K1MlCandidateId,
                  k1ml::K1TruthStatus, k1ml::K1TruthChannel,
                  k1ml::K1MotherPdg, k1ml::K1MotherId,
                  k1ml::K1GeneratedPt, k1ml::K1GeneratedY);
DECLARE_SOA_TABLE(K1MlGenAudit, "AOD", "K1MLGENAUDIT",
                  o2::soa::Index<>, k1ml::K1RecoCollisionId, k1ml::K1OriginalMcParticleId,
                  k1ml::K1MotherPdg, k1ml::K1DaughterPdg1, k1ml::K1DaughterPdg2,
                  k1ml::K1TruthChannel, k1ml::K1GeneratedPt, k1ml::K1GeneratedY,
                  k1ml::K1GenSelectedRecoEvent);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFK1MLTABLES_H_
