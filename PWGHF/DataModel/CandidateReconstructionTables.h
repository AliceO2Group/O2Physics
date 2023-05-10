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

/// \file CandidateReconstructionTables.h
/// \brief Definitions of tables produced by candidate reconstruction workflows
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_
#define PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_

#include <Math/Vector4D.h>
#include <Math/GenVector/Boost.h>

#include "ALICE3/DataModel/ECAL.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::analysis;

namespace o2::aod
{
namespace hf_sel_collision
{
DECLARE_SOA_COLUMN(WhyRejectColl, whyRejectColl, int); //!
} // namespace hf_sel_collision

DECLARE_SOA_TABLE(HfSelCollision, "AOD", "HFSELCOLLISION", //!
                  hf_sel_collision::WhyRejectColl);

namespace hf_sel_track
{
DECLARE_SOA_COLUMN(IsSelProng, isSelProng, int); //!
} // namespace hf_sel_track

DECLARE_SOA_TABLE(HfSelTrack, "AOD", "HFSELTRACK", //!
                  hf_sel_track::IsSelProng);

namespace hf_pv_refit_track
{
DECLARE_SOA_COLUMN(PvRefitX, pvRefitX, float);             //!
DECLARE_SOA_COLUMN(PvRefitY, pvRefitY, float);             //!
DECLARE_SOA_COLUMN(PvRefitZ, pvRefitZ, float);             //!
DECLARE_SOA_COLUMN(PvRefitSigmaX2, pvRefitSigmaX2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXY, pvRefitSigmaXY, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaY2, pvRefitSigmaY2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXZ, pvRefitSigmaXZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaYZ, pvRefitSigmaYZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaZ2, pvRefitSigmaZ2, float); //!
DECLARE_SOA_COLUMN(PvRefitDcaXY, pvRefitDcaXY, float);     //!
DECLARE_SOA_COLUMN(PvRefitDcaZ, pvRefitDcaZ, float);       //!
} // namespace hf_pv_refit_track

DECLARE_SOA_TABLE(HfPvRefitTrack, "AOD", "HFPVREFITTRACK", //!
                  hf_pv_refit_track::PvRefitX,
                  hf_pv_refit_track::PvRefitY,
                  hf_pv_refit_track::PvRefitZ,
                  hf_pv_refit_track::PvRefitSigmaX2,
                  hf_pv_refit_track::PvRefitSigmaXY,
                  hf_pv_refit_track::PvRefitSigmaY2,
                  hf_pv_refit_track::PvRefitSigmaXZ,
                  hf_pv_refit_track::PvRefitSigmaYZ,
                  hf_pv_refit_track::PvRefitSigmaZ2,
                  hf_pv_refit_track::PvRefitDcaXY,
                  hf_pv_refit_track::PvRefitDcaZ);

using BigTracks = soa::Join<Tracks, TracksCov, TracksExtra>;
using BigTracksExtended = soa::Join<BigTracks, aod::TracksDCA>;
using BigTracksMC = soa::Join<BigTracks, McTrackLabels>;
using BigTracksPID = soa::Join<BigTracks,
                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using BigTracksPIDExtended = soa::Join<BigTracksPID, aod::TracksDCA>;

// FIXME: this is a workaround until we get the index columns to work with joins.

namespace hf_track_index
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, Tracks, "_0"); //! Index to first prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, Tracks, "_1"); //! Index to second prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, Tracks, "_2"); //! Index to third prong
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                 //! Index to V0 prong
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);                       //! Index to cascade prong
DECLARE_SOA_COLUMN(HFflag, hfflag, uint8_t);                      //!

DECLARE_SOA_COLUMN(FlagD0ToKPi, flagD0ToKPi, uint8_t);       //!
DECLARE_SOA_COLUMN(FlagJpsiToEE, flagJpsiToEE, uint8_t);     //!
DECLARE_SOA_COLUMN(FlagJpsiToMuMu, flagJpsiToMuMu, uint8_t); //!

DECLARE_SOA_COLUMN(FlagDplusToPiKPi, flagDplusToPiKPi, uint8_t); //!
DECLARE_SOA_COLUMN(FlagLcToPKPi, flagLcToPKPi, uint8_t);         //!
DECLARE_SOA_COLUMN(FlagDsToKKPi, flagDsToKKPi, uint8_t);         //!
DECLARE_SOA_COLUMN(FlagXicToPKPi, flagXicToPKPi, uint8_t);       //!
} // namespace hf_track_index

DECLARE_SOA_TABLE(Hf2Prongs_000, "AOD", "HF2PRONG", //! Table for HF 2 prong candidates (Run 2 converted format)
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_track_index::HFflag);

DECLARE_SOA_TABLE_VERSIONED(Hf2Prongs_001, "AOD", "HF2PRONG", 1, //! Table for HF 2 prong candidates (Run 3 format)
                            o2::soa::Index<>,
                            hf_track_index::CollisionId,
                            hf_track_index::Prong0Id,
                            hf_track_index::Prong1Id,
                            hf_track_index::HFflag);

using Hf2Prongs = Hf2Prongs_001;
using Hf2Prong = Hf2Prongs::iterator;

DECLARE_SOA_TABLE(HfCascades_000, "AOD", "HFCASCADE", //! Table for HF candidates with a V0 (Run 2 converted format)
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::V0Id);

DECLARE_SOA_TABLE_VERSIONED(HfCascades_001, "AOD", "HFCASCADE", 1, //! Table for HF candidates with a V0 (Run 3 format)
                            o2::soa::Index<>,
                            hf_track_index::CollisionId,
                            hf_track_index::Prong0Id,
                            hf_track_index::V0Id);

using HfCascades = HfCascades_001;
using HfCascade = HfCascades::iterator;

DECLARE_SOA_TABLE(Hf3Prongs_000, "AOD", "HF3PRONG", //! Table for HF 3 prong candidates (Run 2 converted format)
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_track_index::Prong2Id,
                  hf_track_index::HFflag);

DECLARE_SOA_TABLE_VERSIONED(Hf3Prongs_001, "AOD", "HF3PRONG", 1, //! Table for HF 3 prong candidates (Run 3 format)
                            o2::soa::Index<>,
                            hf_track_index::CollisionId,
                            hf_track_index::Prong0Id,
                            hf_track_index::Prong1Id,
                            hf_track_index::Prong2Id,
                            hf_track_index::HFflag);

using Hf3Prongs = Hf3Prongs_001;
using Hf3Prong = Hf3Prongs::iterator;

DECLARE_SOA_TABLE(HfCascLf2Prongs, "AOD", "HFCASCLF2PRONG", //! Table for HF 2 prong candidates with a Cascade
                  o2::soa::Index<>,
                  hf_track_index::CascadeId,
                  hf_track_index::Prong0Id);
using HfCascLf2Prong = HfCascLf2Prongs::iterator;

DECLARE_SOA_TABLE(HfCascLf3Prongs, "AOD", "HFCASCLF3PRONG", //! Table for HF 3 prong candidates with a Cascade
                  o2::soa::Index<>,
                  hf_track_index::CascadeId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);
using HfCascLf3Prong = HfCascLf3Prongs::iterator;

namespace hf_track_index
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngD0, prongD0, int, Hf2Prongs, ""); //! Index to a D0 prong
} // namespace hf_track_index

DECLARE_SOA_TABLE(HfDstars, "AOD", "HFDSTAR", //! D* -> D0pi candidates
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::ProngD0Id);
using HfDstar = HfDstars::iterator;

DECLARE_SOA_TABLE(HfCutStatus2Prong, "AOD", "HFCUTSTATUS2P", //!
                  hf_track_index::FlagD0ToKPi,
                  hf_track_index::FlagJpsiToEE,
                  hf_track_index::FlagJpsiToMuMu);

DECLARE_SOA_TABLE(HfCutStatus3Prong, "AOD", "HFCUTSTATUS3P", //!
                  hf_track_index::FlagDplusToPiKPi,
                  hf_track_index::FlagLcToPKPi,
                  hf_track_index::FlagDsToKKPi,
                  hf_track_index::FlagXicToPKPi);

namespace hf_pv_refit_cand_2prong
{
DECLARE_SOA_COLUMN(PvRefitX, pvRefitX, float);             //!
DECLARE_SOA_COLUMN(PvRefitY, pvRefitY, float);             //!
DECLARE_SOA_COLUMN(PvRefitZ, pvRefitZ, float);             //!
DECLARE_SOA_COLUMN(PvRefitSigmaX2, pvRefitSigmaX2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXY, pvRefitSigmaXY, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaY2, pvRefitSigmaY2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXZ, pvRefitSigmaXZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaYZ, pvRefitSigmaYZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaZ2, pvRefitSigmaZ2, float); //!
} // namespace hf_pv_refit_cand_2prong

DECLARE_SOA_TABLE(HfPvRefit2Prong, "AOD", "HFPVREFIT2PRONG", //!
                  hf_pv_refit_cand_2prong::PvRefitX,
                  hf_pv_refit_cand_2prong::PvRefitY,
                  hf_pv_refit_cand_2prong::PvRefitZ,
                  hf_pv_refit_cand_2prong::PvRefitSigmaX2,
                  hf_pv_refit_cand_2prong::PvRefitSigmaXY,
                  hf_pv_refit_cand_2prong::PvRefitSigmaY2,
                  hf_pv_refit_cand_2prong::PvRefitSigmaXZ,
                  hf_pv_refit_cand_2prong::PvRefitSigmaYZ,
                  hf_pv_refit_cand_2prong::PvRefitSigmaZ2);

namespace hf_pv_refit_cand_3prong
{
DECLARE_SOA_COLUMN(PvRefitX, pvRefitX, float);             //!
DECLARE_SOA_COLUMN(PvRefitY, pvRefitY, float);             //!
DECLARE_SOA_COLUMN(PvRefitZ, pvRefitZ, float);             //!
DECLARE_SOA_COLUMN(PvRefitSigmaX2, pvRefitSigmaX2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXY, pvRefitSigmaXY, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaY2, pvRefitSigmaY2, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaXZ, pvRefitSigmaXZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaYZ, pvRefitSigmaYZ, float); //!
DECLARE_SOA_COLUMN(PvRefitSigmaZ2, pvRefitSigmaZ2, float); //!
} // namespace hf_pv_refit_cand_3prong

DECLARE_SOA_TABLE(HfPvRefit3Prong, "AOD", "HFPVREFIT3PRONG", //!
                  hf_pv_refit_cand_3prong::PvRefitX,
                  hf_pv_refit_cand_3prong::PvRefitY,
                  hf_pv_refit_cand_3prong::PvRefitZ,
                  hf_pv_refit_cand_3prong::PvRefitSigmaX2,
                  hf_pv_refit_cand_3prong::PvRefitSigmaXY,
                  hf_pv_refit_cand_3prong::PvRefitSigmaY2,
                  hf_pv_refit_cand_3prong::PvRefitSigmaXZ,
                  hf_pv_refit_cand_3prong::PvRefitSigmaYZ,
                  hf_pv_refit_cand_3prong::PvRefitSigmaZ2);

// general decay properties
namespace hf_cand
{
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float); //!
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float); //!
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, //!
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float); //! sum of (non-weighted) distances of the secondary vertex to its prongs
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //!
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong0, pt2Prong0, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng0, pVectorProng0, //!
                           [](float px, float py, float pz) -> array<float, 3> { return array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter0, errorImpactParameter0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong1, pt2Prong1, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng1, pVectorProng1, //!
                           [](float px, float py, float pz) -> array<float, 3> { return array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter1, errorImpactParameter1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float); //!
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float); //!
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng2, ptProng2, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong2, pt2Prong2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng2, pVectorProng2, //!
                           [](float px, float py, float pz) -> array<float, 3> { return array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter2, impactParameter2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter2, errorImpactParameter2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2, pt2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(P2, p2, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p2(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, //!
                           [](float px, float py, float pz) -> array<float, 3> { return array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //!
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::y(array{px, py, pz}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E2, e2, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e2(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);     //!
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa,                               //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}, array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CPAXY, cpaXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}, array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(Ct, ct, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(array{px, py, pz}, RecoDecay::distance(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}, array{px, py, pz}); });
} // namespace hf_cand

// specific 2-prong decay properties
namespace hf_cand_2prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, //!
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m2(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStar, cosThetaStar, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m, double mTot, int iProng) -> float { return RecoDecay::cosThetaStar(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m, mTot, iProng); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}, errDlxy, array{pxM, pyM}, array{ip0, ip1}, array{errIp0, errIp1}, array{array{px0, py0}, array{px1, py1}}); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level

// mapping of decay types
enum DecayType { D0ToPiK = 0,
                 JpsiToEE,
                 JpsiToMuMu,
                 N2ProngDecays }; // always keep N2ProngDecays at the end

// functions for specific particles

// D0(bar) → π± K∓

template <typename T>
auto ctD0(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kD0));
}

template <typename T>
auto yD0(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kD0));
}

template <typename T>
auto eD0(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kD0));
}

template <typename T>
auto invMassD0ToPiK(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus)});
}

template <typename T>
auto invMassD0barToKPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto cosThetaStarD0(const T& candidate)
{
  return candidate.cosThetaStar(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus)}, RecoDecay::getMassPDG(pdg::Code::kD0), 1);
}

template <typename T>
auto cosThetaStarD0bar(const T& candidate)
{
  return candidate.cosThetaStar(array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)}, RecoDecay::getMassPDG(pdg::Code::kD0), 0);
}

// J/ψ

template <typename T>
auto ctJpsi(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kJPsi));
}

template <typename T>
auto yJpsi(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kJPsi));
}

template <typename T>
auto eJpsi(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kJPsi));
}

// J/ψ → e+ e−
template <typename T>
auto invMassJpsiToEE(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kElectron), RecoDecay::getMassPDG(kElectron)});
}
// J/ψ → μ+ μ−

template <typename T>
auto invMassJpsiToMuMu(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kMuonPlus), RecoDecay::getMassPDG(kMuonMinus)});
}

} // namespace hf_cand_2prong

// general columns
#define HFCAND_COLUMNS                                                                                                                                                                             \
  hf_cand::CollisionId,                                                                                                                                                                            \
    collision::PosX, collision::PosY, collision::PosZ,                                                                                                                                             \
    hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,                                                                                                               \
    hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,                                                                                                                                        \
    hf_cand::Chi2PCA,                                                                                                                                                                              \
    /* dynamic columns */ hf_cand::RSecondaryVertex<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,                                                                                         \
    hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,                                      \
    hf_cand::DecayLengthXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,                                                                                \
    hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ErrorDecayLength>, \
    hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY>,                                         \
    /* prong 0 */ hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,                                                                                  \
    hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,                                                                                                                                       \
    hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,                                                                                                                                      \
    hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,                                                                                                               \
    /* prong 1 */ hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,                                                                                  \
    hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,                                                                                                                                       \
    hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,                                                                                                                                      \
    hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>

// 2-prong decay candidate table
DECLARE_SOA_TABLE(HfCand2ProngBase, "AOD", "HFCAND2PBASE", //!
                  o2::soa::Index<>,
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCand2ProngExt, HfCand2ProngBase, "HFCAND2PEXT", //!
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCand2Prong = HfCand2ProngExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcRec, "AOD", "HFCAND2PMCREC", //!
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcGen, "AOD", "HFCAND2PMCGEN", //!
                  hf_cand_2prong::FlagMcMatchGen,
                  hf_cand_2prong::OriginMcGen);

// cascade decay candidate table

namespace hf_cand_casc
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
// DECLARE_SOA_DYNAMIC_COLUMN(M, m, [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) { return RecoDecay::M(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Pos, ptV0Pos, //!
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Neg, ptV0Neg, //!
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(V0X, v0x, float);
DECLARE_SOA_COLUMN(V0Y, v0y, float);
DECLARE_SOA_COLUMN(V0Z, v0z, float);

template <typename T>
auto invMassLcToK0sP(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kProton), RecoDecay::getMassPDG(kK0Short)}); // first daughter is bachelor
}

template <typename T>
auto invMassGammaToEE(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kElectron), RecoDecay::getMassPDG(kElectron)});
}

} // namespace hf_cand_casc

DECLARE_SOA_TABLE(HfCandCascBase, "AOD", "HFCANDCASCBASE", //!
                                                           // general columns
                  HFCAND_COLUMNS,
                  // cascade specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_track_index::Prong0Id,
                  hf_track_index::V0Id, // V0 index
                  // V0
                  hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z,
                  v0data::PosTrackId, v0data::NegTrackId, // indices of V0 tracks in FullTracks table
                  v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                  v0data::DCAV0Daughters,
                  v0data::DCAPosToPV, // this is the impact param wrt prim vtx in xy!
                  v0data::DCANegToPV, // this is the impact param wrt prim vtx in xy!
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Pt2<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::P<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::P2<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::PVector<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_casc::Px, hf_cand_casc::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::Phi<hf_cand_casc::Px, hf_cand_casc::Py>,
                  hf_cand::Y<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::E<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  hf_cand::E2<hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz>,
                  // dynamic columns from V0
                  hf_cand_casc::PtV0Pos<v0data::PxPos, v0data::PyPos>, // pT of positive V0 daughter
                  hf_cand_casc::PtV0Neg<v0data::PxNeg, v0data::PyNeg>, // pT of negative V0 daughter
                  v0data::V0Radius<hf_cand_casc::V0X, hf_cand_casc::V0Y>,
                  v0data::V0CosPA<hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, collision::PosX, collision::PosY, collision::PosZ>,
                  v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>);
//                  ,
//                  v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
//                  v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
//                  v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandCascExt, HfCandCascBase, "HFCANDCASCEXT", //!
                                hf_cand_casc::Px, hf_cand_casc::Py, hf_cand_casc::Pz);

using HfCandCascade = HfCandCascExt;

// table with results of reconstruction level MC matching for Cascade
DECLARE_SOA_TABLE(HfCandCascadeMcRec, "AOD", "HFCANDCASCMCREC", //!
                  hf_cand_casc::FlagMcMatchRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandCascadeMcGen, "AOD", "HFCANDCASCMCGEN", //!
                  hf_cand_casc::FlagMcMatchGen);

// specific BPlus candidate properties
namespace hf_cand_bplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // D0 index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level

enum DecayType { BplusToD0Pi = 0 };

// B± → D0bar(D0) π±

template <typename T>
auto ctBplus(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kBPlus));
}

template <typename T>
auto yBplus(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kBPlus));
}

template <typename T>
auto eBplus(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kBPlus));
}

template <typename T>
auto invMassBplusToD0Pi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(pdg::Code::kD0), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto cosThetaStarBplus(const T& candidate)
{
  return candidate.cosThetaStar(array{RecoDecay::getMassPDG(pdg::Code::kD0), RecoDecay::getMassPDG(kPiPlus)}, RecoDecay::getMassPDG(pdg::Code::kBPlus), 1);
}
} // namespace hf_cand_bplus

// declare dedicated BPlus decay candidate table
DECLARE_SOA_TABLE(HfCandBplusBase, "AOD", "HFCANDBPLUSBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  o2::soa::Index<>,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_bplus::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBplusExt, HfCandBplusBase, "HFCANDBPLUSEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandBplus = HfCandBplusExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcRec, "AOD", "HFCANDBPMCREC",
                  hf_cand_bplus::FlagMcMatchRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcGen, "AOD", "HFCANDBPMCGEN",
                  hf_cand_bplus::FlagMcMatchGen);

// specific 3-prong decay properties
namespace hf_cand_3prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1 + 1.f * aod::hf_cand::pxProng2);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1 + 1.f * aod::hf_cand::pyProng2);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1 + 1.f * aod::hf_cand::pzProng2);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const array<double, 3>& m) -> float { return RecoDecay::m(array{array{px0, py0, pz0}, array{px1, py1, pz1}, array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const array<double, 3>& m) -> float { return RecoDecay::m2(array{array{px0, py0, pz0}, array{px1, py1, pz1}, array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1, float impParProng2) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1, impParProng2); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float ip2, float errIp2, float px0, float py0, float px1, float py1, float px2, float py2) -> float { return RecoDecay::maxNormalisedDeltaIP(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}, errDlxy, array{pxM, pyM}, array{ip0, ip1, ip2}, array{errIp0, errIp1, errIp2}, array{array{px0, py0}, array{px1, py1}, array{px2, py2}}); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t); //! swapping of the prongs order
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, generator level

// mapping of decay types
enum DecayType { DplusToPiKPi = 0,
                 LcToPKPi,
                 DsToKKPi,
                 XicToPKPi,
                 N3ProngDecays }; // always keep N3ProngDecays at the end

// functions for specific particles

// D± → π± K∓ π±

template <typename T>
auto ctDplus(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kDPlus));
}

template <typename T>
auto yDplus(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kDPlus));
}

template <typename T>
auto eDplus(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kDPlus));
}

template <typename T>
auto invMassDplusToPiKPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)});
}

// Ds± → K± K∓ π±

template <typename T>
auto ctDs(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kDS));
}

template <typename T>
auto yDs(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kDS));
}

template <typename T>
auto eDs(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kDS));
}

template <typename T>
auto invMassDsToKKPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto invMassDsToPiKK(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kKPlus)});
}

template <typename T>
auto deltaMassPhiDsToKKPi(const T& candidate)
{
  double invMassKKpair = RecoDecay::m(array{candidate.pVectorProng0(), candidate.pVectorProng1()}, array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kKPlus)});
  return std::abs(invMassKKpair - RecoDecay::getMassPDG(pdg::Code::kPhi));
}

template <typename T>
auto deltaMassPhiDsToPiKK(const T& candidate)
{
  double invMassKKpair = RecoDecay::m(array{candidate.pVectorProng1(), candidate.pVectorProng2()}, array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kKPlus)});
  return std::abs(invMassKKpair - RecoDecay::getMassPDG(pdg::Code::kPhi));
}

/// Calculate the cosine of the angle between the pion and the opposite sign kaon in the phi rest frame
/// \param candidate Ds candidate from aod::HfCand3Prong table
/// \param option mass hypothesis considered: 0 = KKPi, 1 = PiKK
/// \return cosine of pion-kaon angle in the phi rest frame
template <typename T>
auto cosPiKPhiRestFrame(const T& candidate, int option)
{
  // Ported from AliAODRecoDecayHF3Prong::CosPiKPhiRFrame
  array<float, 3> momPi;
  array<float, 3> momK1;
  array<float, 3> momK2;

  if (option == 0) { // KKPi
    momPi = candidate.pVectorProng2();
    momK1 = candidate.pVectorProng1();
    momK2 = candidate.pVectorProng0();
  } else { // PiKK
    momPi = candidate.pVectorProng0();
    momK1 = candidate.pVectorProng1();
    momK2 = candidate.pVectorProng2();
  }

  ROOT::Math::PxPyPzMVector vecPi(momPi[0], momPi[1], momPi[2], RecoDecay::getMassPDG(kPiPlus));
  ROOT::Math::PxPyPzMVector vecK1(momK1[0], momK1[1], momK1[2], RecoDecay::getMassPDG(kKPlus));
  ROOT::Math::PxPyPzMVector vecK2(momK2[0], momK2[1], momK2[2], RecoDecay::getMassPDG(kKPlus));
  ROOT::Math::PxPyPzMVector vecPhi = vecK1 + vecK2;

  ROOT::Math::Boost boostToPhiRestFrame(vecPhi.BoostToCM());
  auto momPiPhiRestFrame = boostToPhiRestFrame(vecPi).Vect();
  auto momK1PhiRestFrame = boostToPhiRestFrame(vecK1).Vect();

  return momPiPhiRestFrame.Dot(momK1PhiRestFrame) / std::sqrt(momPiPhiRestFrame.Mag2() * momK1PhiRestFrame.Mag2());
}

template <typename T>
auto cos3PiKDsToKKPi(const T& candidate)
{
  auto cosPiK = cosPiKPhiRestFrame(candidate, 0);
  return cosPiK * cosPiK * cosPiK;
}

template <typename T>
auto cos3PiKDsToPiKK(const T& candidate)
{
  auto cosPiK = cosPiKPhiRestFrame(candidate, 1);
  return cosPiK * cosPiK * cosPiK;
}

// Λc± → p± K∓ π±

template <typename T>
auto ctLc(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus));
}

template <typename T>
auto yLc(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus));
}

template <typename T>
auto eLc(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus));
}

template <typename T>
auto invMassLcToPKPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kProton), RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto invMassLcToPiKP(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kProton)});
}

// Ξc± → p± K∓ π±

template <typename T>
auto ctXic(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kXiCPlus));
}

template <typename T>
auto yXic(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kXiCPlus));
}

template <typename T>
auto eXic(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kXiCPlus));
}

template <typename T>
auto invMassXicToPKPi(const T& candidate)
{
  return invMassLcToPKPi(candidate);
}

template <typename T>
auto invMassXicToPiKP(const T& candidate)
{
  return invMassLcToPiKP(candidate);
}

} // namespace hf_cand_3prong

// 3-prong decay candidate table
DECLARE_SOA_TABLE(HfCand3ProngBase, "AOD", "HFCAND3PBASE", //!
                  o2::soa::Index<>,
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2>,
                  /* prong 2 */
                  hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCand3ProngExt, HfCand3ProngBase, "HFCAND3PEXT", //!
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

using HfCand3Prong = HfCand3ProngExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcRec, "AOD", "HFCAND3PMCREC", //!
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::IsCandidateSwapped,
                  hf_cand_3prong::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcGen, "AOD", "HFCAND3PMCGEN", //!
                  hf_cand_3prong::FlagMcMatchGen,
                  hf_cand_3prong::OriginMcGen,
                  hf_cand_3prong::FlagMcDecayChanGen);

namespace hf_cand_casc_lf_2prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, //!
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m2(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStar, cosThetaStar, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m, double mTot, int iProng) -> float { return RecoDecay::cosThetaStar(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m, mTot, iProng); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}, errDlxy, array{pxM, pyM}, array{ip0, ip1}, array{errIp0, errIp1}, array{array{px0, py0}, array{px1, py1}}); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level

template <typename T>
auto invMassXiczeroToXiPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kXiMinus), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto invMassOmegaczeroToOmegaPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kOmegaMinus), RecoDecay::getMassPDG(kPiPlus)});
}

// mapping of decay types
enum DecayType { XiczeroToXiPi = 0,
                 OmegaczeroToOmegaPi,
                 N2ProngDecays }; // always keep N2ProngDecays at the end

} // namespace hf_cand_casc_lf_2prong

namespace hf_cand_casc_lf_3prong
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProduct, impactParameterProduct, //!
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::m2(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStar, cosThetaStar, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m, double mTot, int iProng) -> float { return RecoDecay::cosThetaStar(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m, mTot, iProng); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(array{xVtxP, yVtxP}, array{xVtxS, yVtxS}, errDlxy, array{pxM, pyM}, array{ip0, ip1}, array{errIp0, errIp1}, array{array{px0, py0}, array{px1, py1}}); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level

template <typename T>
auto invMassXicplusToXiPiPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kXiMinus), RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kPiPlus)});
}

// mapping of decay types
enum DecayType { XicplusToXiPiPi = 0,
                 N3ProngDecays }; // always keep N3ProngDecays at the end

} // namespace hf_cand_casc_lf_3prong

namespace hf_cand_x
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // Jpsi index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel flag, generator level

// mapping of decay types
enum DecayType { XToJpsiToEEPiPi = 0,
                 XToJpsiToMuMuPiPi }; // move this to a dedicated cascade namespace in the future?

// X → Jpsi π+ π-
// TODO: add pdg code for X (9920443), temporarily hardcode mass here:
float massX = 3.872; // replace this with: "RecoDecay::getMassPDG(9920443)" when pdg is added
template <typename T>
auto ctX(const T& candidate)
{
  return candidate.ct(massX);
}

template <typename T>
auto yX(const T& candidate)
{
  return candidate.y(massX);
}

template <typename T>
auto eX(const T& candidate)
{
  return candidate.e(massX);
}

template <typename T>
auto invMassXToJpsiPiPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(443), RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kPiPlus)});
}

/// Difference between the X mass and the sum of the J/psi and di-pion masses
template <typename T>
auto qX(const T& candidate)
{
  auto piVec1 = array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()};
  auto piVec2 = array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()};
  double massPi = RecoDecay::getMassPDG(kPiPlus);

  auto arrayMomenta = array{piVec1, piVec2};
  double massPiPi = RecoDecay::m(arrayMomenta, array{massPi, massPi});

  // PDG mass, as reported in CMS paper https://arxiv.org/pdf/1302.3968.pdf
  double massJpsi = RecoDecay::getMassPDG(o2::analysis::pdg::kJPsi);

  double massX = invMassXToJpsiPiPi(candidate);
  return std::abs(massX - massJpsi - massPiPi);
}

/// Angular difference between the J/psi and the pion
template <typename T>
auto dRX(const T& candidate, int numPi)
{
  double etaJpsi = RecoDecay::eta(array{candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()});
  double phiJpsi = RecoDecay::phi(candidate.pxProng0(), candidate.pyProng0());

  double etaPi, phiPi;

  if (numPi <= 1) {
    etaPi = RecoDecay::eta(array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()});
    phiPi = RecoDecay::phi(candidate.pxProng1(), candidate.pyProng1());
  } else {
    etaPi = RecoDecay::eta(array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()});
    phiPi = RecoDecay::phi(candidate.pxProng2(), candidate.pyProng2());
  }

  double deltaEta = etaJpsi - etaPi;
  double deltaPhi = RecoDecay::constrainAngle(phiJpsi - phiPi, -o2::constants::math::PI);

  return RecoDecay::sqrtSumOfSquares(deltaEta, deltaPhi);
}

/// Difference in pT between the two pions
template <typename T>
auto balancePtPionsX(const T& candidate)
{
  double ptPi1 = RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1());
  double ptPi2 = RecoDecay::pt(candidate.pxProng2(), candidate.pyProng2());
  return std::abs(ptPi1 - ptPi2) / (ptPi1 + ptPi2);
}
} // namespace hf_cand_x

// declare dedicated X candidate table
DECLARE_SOA_TABLE(HfCandXBase, "AOD", "HFCANDXBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  hf_cand_x::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id, // note the difference between Jpsi and pion indices
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_3prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand_3prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1, hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* prong 2 */
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Pt2Prong2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::PVectorProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Pt2<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::P2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::E2<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXExt, HfCandXBase, "HFCANDXEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

using HfCandX = HfCandXExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXMcRec, "AOD", "HFCANDXMCREC", //!
                  hf_cand_x::FlagMcMatchRec,
                  hf_cand_x::OriginMcRec,
                  hf_cand_x::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXMcGen, "AOD", "HFCANDXMCGEN", //!
                  hf_cand_x::FlagMcMatchGen,
                  hf_cand_x::OriginMcGen,
                  hf_cand_x::FlagMcDecayChanGen);

// definition of columns and tables for D-Dbar correlation pairs
namespace hf_correlation_d_dbar
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtD, ptD, float);
DECLARE_SOA_COLUMN(PtDbar, ptDbar, float);
DECLARE_SOA_COLUMN(MD, mD, float);
DECLARE_SOA_COLUMN(MDbar, mDbar, float);
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int);
} // namespace hf_correlation_d_dbar
DECLARE_SOA_TABLE(DDbarPair, "AOD", "DDBARPAIR",
                  aod::hf_correlation_d_dbar::DeltaPhi,
                  aod::hf_correlation_d_dbar::DeltaEta,
                  aod::hf_correlation_d_dbar::PtD,
                  aod::hf_correlation_d_dbar::PtDbar);
DECLARE_SOA_TABLE(DDbarRecoInfo, "AOD", "DDBARRECOINFO",
                  aod::hf_correlation_d_dbar::MD,
                  aod::hf_correlation_d_dbar::MDbar,
                  aod::hf_correlation_d_dbar::SignalStatus);

// definition of columns and tables for D0-Hadron correlation pairs
namespace hf_correlation_d0_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtD, ptD, float);
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);
DECLARE_SOA_COLUMN(MD, mD, float);
DECLARE_SOA_COLUMN(MDbar, mDbar, float);
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int);
} // namespace hf_correlation_d0_hadron
DECLARE_SOA_TABLE(DHadronPair, "AOD", "DHADRONPAIR",
                  aod::hf_correlation_d0_hadron::DeltaPhi,
                  aod::hf_correlation_d0_hadron::DeltaEta,
                  aod::hf_correlation_d0_hadron::PtD,
                  aod::hf_correlation_d0_hadron::PtHadron);
DECLARE_SOA_TABLE(DHadronRecoInfo, "AOD", "DHADRONRECOINFO",
                  aod::hf_correlation_d0_hadron::MD,
                  aod::hf_correlation_d0_hadron::MDbar,
                  aod::hf_correlation_d0_hadron::SignalStatus);

// definition of columns and tables for Ds-Hadron correlation pairs
namespace hf_correlation_ds_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);        //! DeltaPhi between Ds and Hadrons
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);        //! DeltaEta between Ds and Hadrons
DECLARE_SOA_COLUMN(PtD, ptD, float);                  //! Transverse momentum of Ds
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);        //! Transverse momentum of Hadron
DECLARE_SOA_COLUMN(MD, mD, float);                    //! Invariant mass of Ds
DECLARE_SOA_COLUMN(PoolBin, poolBin, int);            //! Pool Bin for the MixedEvent
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, bool); //! Used in MC-Rec, Ds Signal
DECLARE_SOA_COLUMN(PromptStatus, promptStatus, bool); //! Used in MC-Rec, Ds Prompt
} // namespace hf_correlation_ds_hadron
DECLARE_SOA_TABLE(DsHadronPair, "AOD", "DSHPAIR", //! Ds-Hadrons pairs Informations
                  aod::hf_correlation_ds_hadron::DeltaPhi,
                  aod::hf_correlation_ds_hadron::DeltaEta,
                  aod::hf_correlation_ds_hadron::PtD,
                  aod::hf_correlation_ds_hadron::PtHadron,
                  aod::hf_correlation_ds_hadron::PoolBin);
DECLARE_SOA_TABLE(DsHadronRecoInfo, "AOD", "DSHRECOINFO", //! Ds-Hadrons pairs Reconstructed Informations
                  aod::hf_correlation_ds_hadron::MD,
                  aod::hf_correlation_ds_hadron::SignalStatus,
                  aod::hf_correlation_ds_hadron::PromptStatus);
DECLARE_SOA_TABLE(DsHadronGenInfo, "AOD", "DSHGENINFO", //! Ds-Hadrons pairs Generated Informations
                  aod::hf_correlation_ds_hadron::PromptStatus);

// table for selection of collisions with at least one Ds meson
namespace hf_sel_collision_ds
{
DECLARE_SOA_COLUMN(DsFound, dsFound, bool); //! Ds found in a collision
} // namespace hf_sel_collision_ds
DECLARE_SOA_TABLE(DsSelCollision, "AOD", "DSCOLL", aod::hf_sel_collision_ds::DsFound);

// definition of columns and tables for Dplus-Hadron correlation pairs
namespace hf_correlation_dplus_hadron
{
DECLARE_SOA_COLUMN(DeltaPhi, deltaPhi, float);
DECLARE_SOA_COLUMN(DeltaEta, deltaEta, float);
DECLARE_SOA_COLUMN(PtD, ptD, float);
DECLARE_SOA_COLUMN(PtHadron, ptHadron, float);
DECLARE_SOA_COLUMN(MD, mD, float);
DECLARE_SOA_COLUMN(SignalStatus, signalStatus, int);
} // namespace hf_correlation_dplus_hadron
DECLARE_SOA_TABLE(DplusHadronPair, "AOD", "DPLUSHPAIR",
                  aod::hf_correlation_dplus_hadron::DeltaPhi,
                  aod::hf_correlation_dplus_hadron::DeltaEta,
                  aod::hf_correlation_dplus_hadron::PtD,
                  aod::hf_correlation_dplus_hadron::PtHadron);
DECLARE_SOA_TABLE(DplusHadronRecoInfo, "AOD", "DPLUSHRECOINFO",
                  aod::hf_correlation_dplus_hadron::MD,
                  aod::hf_correlation_dplus_hadron::SignalStatus);

// specific Xicc candidate properties
namespace hf_cand_xicc
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Xic index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
// mapping of decay types
enum DecayType { XiccToXicPi = 0 }; // move this to a dedicated cascade namespace in the future?

// Ξcc±± → p± K∓ π± π±

template <typename T>
auto ctXicc(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kXiCCPlusPlus));
}

template <typename T>
auto yXicc(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kXiCCPlusPlus));
}

template <typename T>
auto eXicc(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kXiCCPlusPlus));
}

template <typename T>
auto invMassXiccToXicPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(pdg::Code::kXiCPlus), RecoDecay::getMassPDG(kPiPlus)});
}
} // namespace hf_cand_xicc

// declare dedicated Xicc candidate table
DECLARE_SOA_TABLE(HfCandXiccBase, "AOD", "HFCANDXICCBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_xicc::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXiccExt, HfCandXiccBase, "HFCANDXICCEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandXicc = HfCandXiccExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXiccMcRec, "AOD", "HFCANDXICCMCREC", //!
                  hf_cand_xicc::FlagMcMatchRec,
                  hf_cand_xicc::OriginMcRec,
                  hf_cand_xicc::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXiccMcGen, "AOD", "HFCANDXICCMCGEN", //!
                  hf_cand_xicc::FlagMcMatchGen,
                  hf_cand_xicc::OriginMcGen);

// specific Omegac and Xic to Xi Pi candidate properties
namespace hf_cand_toxipi
{
// Data processing results:
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(XPv, xPv, float);
DECLARE_SOA_COLUMN(YPv, yPv, float);
DECLARE_SOA_COLUMN(ZPv, zPv, float);
DECLARE_SOA_COLUMN(XDecayVtxOmegac, xDecayVtxOmegac, float);
DECLARE_SOA_COLUMN(YDecayVtxOmegac, yDecayVtxOmegac, float);
DECLARE_SOA_COLUMN(ZDecayVtxOmegac, zDecayVtxOmegac, float);
DECLARE_SOA_COLUMN(XDecayVtxCascade, xDecayVtxCascade, float);
DECLARE_SOA_COLUMN(YDecayVtxCascade, yDecayVtxCascade, float);
DECLARE_SOA_COLUMN(ZDecayVtxCascade, zDecayVtxCascade, float);
DECLARE_SOA_COLUMN(XDecayVtxV0, xDecayVtxV0, float);
DECLARE_SOA_COLUMN(YDecayVtxV0, yDecayVtxV0, float);
DECLARE_SOA_COLUMN(ZDecayVtxV0, zDecayVtxV0, float);
DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t); // sign of pi <- xi
DECLARE_SOA_COLUMN(Chi2PCAOmegac, chi2PcaOmegac, float);
DECLARE_SOA_COLUMN(CovVtxOmegac0, covVtxOmegac0, float);
DECLARE_SOA_COLUMN(CovVtxOmegac1, covVtxOmegac1, float);
DECLARE_SOA_COLUMN(CovVtxOmegac2, covVtxOmegac2, float);
DECLARE_SOA_COLUMN(CovVtxOmegac3, covVtxOmegac3, float);
DECLARE_SOA_COLUMN(CovVtxOmegac4, covVtxOmegac4, float);
DECLARE_SOA_COLUMN(CovVtxOmegac5, covVtxOmegac5, float);
DECLARE_SOA_COLUMN(CovVtxV00, covVtxV00, float);
DECLARE_SOA_COLUMN(CovVtxV01, covVtxV01, float);
DECLARE_SOA_COLUMN(CovVtxV02, covVtxV02, float);
DECLARE_SOA_COLUMN(CovVtxV03, covVtxV03, float);
DECLARE_SOA_COLUMN(CovVtxV04, covVtxV04, float);
DECLARE_SOA_COLUMN(CovVtxV05, covVtxV05, float);
DECLARE_SOA_COLUMN(CovVtxCasc0, covVtxCasc0, float);
DECLARE_SOA_COLUMN(CovVtxCasc1, covVtxCasc1, float);
DECLARE_SOA_COLUMN(CovVtxCasc2, covVtxCasc2, float);
DECLARE_SOA_COLUMN(CovVtxCasc3, covVtxCasc3, float);
DECLARE_SOA_COLUMN(CovVtxCasc4, covVtxCasc4, float);
DECLARE_SOA_COLUMN(CovVtxCasc5, covVtxCasc5, float);
DECLARE_SOA_COLUMN(PxOmegac, pxOmegac, float);
DECLARE_SOA_COLUMN(PyOmegac, pyOmegac, float);
DECLARE_SOA_COLUMN(PzOmegac, pzOmegac, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PxPrimaryPi, pxPrimaryPi, float);
DECLARE_SOA_COLUMN(PyPrimaryPi, pyPrimaryPi, float);
DECLARE_SOA_COLUMN(PzPrimaryPi, pzPrimaryPi, float);
DECLARE_SOA_COLUMN(PxLambda, pxLambda, float);
DECLARE_SOA_COLUMN(PyLambda, pyLambda, float);
DECLARE_SOA_COLUMN(PzLambda, pzLambda, float);
DECLARE_SOA_COLUMN(PxPiFromCasc, pxPiFromCasc, float);
DECLARE_SOA_COLUMN(PyPiFromCasc, pyPiFromCasc, float);
DECLARE_SOA_COLUMN(PzPiFromCasc, pzPiFromCasc, float);
DECLARE_SOA_COLUMN(PxPosV0Dau, pxPosV0Dau, float);
DECLARE_SOA_COLUMN(PyPosV0Dau, pyPosV0Dau, float);
DECLARE_SOA_COLUMN(PzPosV0Dau, pzPosV0Dau, float);
DECLARE_SOA_COLUMN(PxNegV0Dau, pxNegV0Dau, float);
DECLARE_SOA_COLUMN(PyNegV0Dau, pyNegV0Dau, float);
DECLARE_SOA_COLUMN(PzNegV0Dau, pzNegV0Dau, float);
DECLARE_SOA_COLUMN(ImpactParCascXY, impactParCascXY, float);
DECLARE_SOA_COLUMN(ImpactParPrimaryPiXY, impactParPrimaryPiXY, float);
DECLARE_SOA_COLUMN(ImpactParCascZ, impactParCascZ, float);
DECLARE_SOA_COLUMN(ImpactParPrimaryPiZ, impactParPrimaryPiZ, float);
DECLARE_SOA_COLUMN(ImpactParV0XY, impactParV0XY, float);
DECLARE_SOA_COLUMN(ImpactParV0Z, impactParV0Z, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParPrimaryPiXY, errImpactParPrimaryPiXY, float);
DECLARE_SOA_COLUMN(ErrImpactParV0XY, errImpactParV0XY, float);
DECLARE_SOA_INDEX_COLUMN(V0, v0);
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);
DECLARE_SOA_INDEX_COLUMN_FULL(PrimaryPi, primaryPi, int, Tracks, "_primarypi");
DECLARE_SOA_COLUMN(ImpactParOmegacXY, impactParOmegacXY, float);
DECLARE_SOA_COLUMN(ImpactParOmegacZ, impactParOmegacZ, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, double);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, double);
DECLARE_SOA_COLUMN(InvMassOmegac, invMassOmegac, double);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, double);
DECLARE_SOA_COLUMN(CosPAOmegac, cosPAOmegac, double);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, double);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, double);
DECLARE_SOA_COLUMN(CosPAXYOmegac, cosPAXYOmegac, double);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, double);
DECLARE_SOA_COLUMN(CTauOmegac, ctauOmegac, double);
DECLARE_SOA_COLUMN(CTauCascade, ctauCascade, double);
DECLARE_SOA_COLUMN(CTauV0, ctauV0, double);
DECLARE_SOA_COLUMN(CTauXic, ctauXic, double);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, double);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, double);
DECLARE_SOA_COLUMN(EtaPiFromCasc, etaPiFromCasc, double);
DECLARE_SOA_COLUMN(EtaPiFromOme, etaPiFromOme, double);
DECLARE_SOA_COLUMN(EtaOmegac, etaOmegac, double);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, double);
DECLARE_SOA_COLUMN(EtaV0, etaV0, double);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau0, dcaZToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau1, dcaZToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaZToPvCascDau, dcaZToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaOmegacDau, dcaOmegacDau, float);

// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level

// mapping of decay types
enum DecayType { DecayToXiPi = 0,
                 OmegaczeroToXiPi,
                 XiczeroToXiPi };

} // end of namespace hf_cand_toxipi

// declare dedicated Omegac and Xic to Xi Pi candidate table
DECLARE_SOA_TABLE(HfCandToXiPi, "AOD", "HFCANDTOXIPI",
                  o2::soa::Index<>,
                  hf_cand_toxipi::CollisionId, hf_cand_toxipi::XPv, hf_cand_toxipi::YPv, hf_cand_toxipi::ZPv,
                  hf_cand_toxipi::XDecayVtxOmegac, hf_cand_toxipi::YDecayVtxOmegac, hf_cand_toxipi::ZDecayVtxOmegac,
                  hf_cand_toxipi::XDecayVtxCascade, hf_cand_toxipi::YDecayVtxCascade, hf_cand_toxipi::ZDecayVtxCascade,
                  hf_cand_toxipi::XDecayVtxV0, hf_cand_toxipi::YDecayVtxV0, hf_cand_toxipi::ZDecayVtxV0,
                  hf_cand_toxipi::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_toxipi::Chi2PCAOmegac,
                  hf_cand_toxipi::CovVtxOmegac0, hf_cand_toxipi::CovVtxOmegac1, hf_cand_toxipi::CovVtxOmegac2, hf_cand_toxipi::CovVtxOmegac3, hf_cand_toxipi::CovVtxOmegac4, hf_cand_toxipi::CovVtxOmegac5,
                  hf_cand_toxipi::CovVtxV00, hf_cand_toxipi::CovVtxV01, hf_cand_toxipi::CovVtxV02, hf_cand_toxipi::CovVtxV03, hf_cand_toxipi::CovVtxV04, hf_cand_toxipi::CovVtxV05,
                  hf_cand_toxipi::CovVtxCasc0, hf_cand_toxipi::CovVtxCasc1, hf_cand_toxipi::CovVtxCasc2, hf_cand_toxipi::CovVtxCasc3, hf_cand_toxipi::CovVtxCasc4, hf_cand_toxipi::CovVtxCasc5,
                  hf_cand_toxipi::PxOmegac, hf_cand_toxipi::PyOmegac, hf_cand_toxipi::PzOmegac,
                  hf_cand_toxipi::PxCasc, hf_cand_toxipi::PyCasc, hf_cand_toxipi::PzCasc,
                  hf_cand_toxipi::PxPrimaryPi, hf_cand_toxipi::PyPrimaryPi, hf_cand_toxipi::PzPrimaryPi,
                  hf_cand_toxipi::PxLambda, hf_cand_toxipi::PyLambda, hf_cand_toxipi::PzLambda,
                  hf_cand_toxipi::PxPiFromCasc, hf_cand_toxipi::PyPiFromCasc, hf_cand_toxipi::PzPiFromCasc,
                  hf_cand_toxipi::PxPosV0Dau, hf_cand_toxipi::PyPosV0Dau, hf_cand_toxipi::PzPosV0Dau,
                  hf_cand_toxipi::PxNegV0Dau, hf_cand_toxipi::PyNegV0Dau, hf_cand_toxipi::PzNegV0Dau,
                  hf_cand_toxipi::ImpactParCascXY, hf_cand_toxipi::ImpactParPrimaryPiXY, hf_cand_toxipi::ImpactParCascZ, hf_cand_toxipi::ImpactParPrimaryPiZ,
                  hf_cand_toxipi::ImpactParV0XY, hf_cand_toxipi::ImpactParV0Z,
                  hf_cand_toxipi::ErrImpactParCascXY, hf_cand_toxipi::ErrImpactParPrimaryPiXY, hf_cand_toxipi::ErrImpactParV0XY,
                  hf_cand_toxipi::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_toxipi::CascadeId, hf_cand_toxipi::PrimaryPiId, cascdata::BachelorId,
                  hf_cand_toxipi::ImpactParOmegacXY, hf_cand_toxipi::ImpactParOmegacZ,
                  hf_cand_toxipi::InvMassLambda, hf_cand_toxipi::InvMassCascade, hf_cand_toxipi::InvMassOmegac,
                  hf_cand_toxipi::CosPAV0, hf_cand_toxipi::CosPAOmegac, hf_cand_toxipi::CosPACasc, hf_cand_toxipi::CosPAXYV0, hf_cand_toxipi::CosPAXYOmegac, hf_cand_toxipi::CosPAXYCasc,
                  hf_cand_toxipi::CTauOmegac, hf_cand_toxipi::CTauCascade, hf_cand_toxipi::CTauV0, hf_cand_toxipi::CTauXic,
                  hf_cand_toxipi::EtaV0PosDau, hf_cand_toxipi::EtaV0NegDau, hf_cand_toxipi::EtaPiFromCasc, hf_cand_toxipi::EtaPiFromOme,
                  hf_cand_toxipi::EtaOmegac, hf_cand_toxipi::EtaCascade, hf_cand_toxipi::EtaV0,
                  hf_cand_toxipi::DcaXYToPvV0Dau0, hf_cand_toxipi::DcaXYToPvV0Dau1, hf_cand_toxipi::DcaXYToPvCascDau,
                  hf_cand_toxipi::DcaZToPvV0Dau0, hf_cand_toxipi::DcaZToPvV0Dau1, hf_cand_toxipi::DcaZToPvCascDau,
                  hf_cand_toxipi::DcaCascDau, hf_cand_toxipi::DcaV0Dau, hf_cand_toxipi::DcaOmegacDau, hf_track_index::HFflag);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfToXiPiMCRec, "AOD", "HFTOXIPIMCREC", //!
                  hf_cand_toxipi::FlagMcMatchRec,
                  hf_cand_toxipi::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfToXiPiMCGen, "AOD", "HFTOXIPIMCGEN", //!
                  hf_cand_toxipi::FlagMcMatchGen);

// specific chic candidate properties
namespace hf_cand_chic
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // Jpsi index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, ECALs, "_1");
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); // resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); // resonant decay channel flag, generator level
DECLARE_SOA_COLUMN(JpsiToMuMuMass, jpsiToMuMuMass, float);          // Jpsi mass
// mapping of decay types
enum DecayType { ChicToJpsiToEEGamma = 0,
                 ChicToJpsiToMuMuGamma }; // move this to a dedicated cascade namespace in the future?
// chic → Jpsi gamma
template <typename T>
auto ctChic(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kChiC1));
}

template <typename T>
auto yChic(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kChiC1));
}

template <typename T>
auto eChic(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kChiC1));
}
template <typename T>
auto invMassChicToJpsiGamma(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(pdg::Code::kJPsi), 0.});
}

} // namespace hf_cand_chic

// declare dedicated chi_c candidate table
DECLARE_SOA_TABLE(HfCandChicBase, "AOD", "HFCANDCHICBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_chic::Prong0Id, hf_cand_chic::Prong1Id,
                  hf_track_index::HFflag, hf_cand_chic::JpsiToMuMuMass,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* prong 2 */
                  //                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  //                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  //                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandChicExt, HfCandChicBase, "HFCANDCHICEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandChic = HfCandChicExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandChicMcRec, "AOD", "HFCANDCHICMCREC", //!
                  hf_cand_chic::FlagMcMatchRec,
                  hf_cand_chic::OriginMcRec,
                  hf_cand_chic::FlagMcDecayChanRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandChicMcGen, "AOD", "HFCANDCHICMCGEN", //!
                  hf_cand_chic::FlagMcMatchGen,
                  hf_cand_chic::OriginMcGen,
                  hf_cand_chic::FlagMcDecayChanGen);

// specific Lb candidate properties
namespace hf_cand_lb
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Lb index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
// mapping of decay types
enum DecayType { LbToLcPi }; // move this to a dedicated cascade namespace in the future?

// Λb → Λc+ π- → p K- π+ π-
// float massLb = RecoDecay::getMassPDG(pdg::Code::kLambdaB0);
template <typename T>
auto ctLb(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kLambdaB0));
}

template <typename T>
auto yLb(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kLambdaB0));
}

template <typename T>
auto eLb(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kLambdaB0));
}
template <typename T>
auto invMassLbToLcPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus), RecoDecay::getMassPDG(kPiPlus)});
}
} // namespace hf_cand_lb

// declare dedicated Lb candidate table
DECLARE_SOA_TABLE(HfCandLbBase, "AOD", "HFCANDLBBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 3-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_lb::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandLbExt, HfCandLbBase, "HFCANDLBEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandLb = HfCandLbExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandLbMcRec, "AOD", "HFCANDLBMCREC", //!
                  hf_cand_lb::FlagMcMatchRec,
                  hf_cand_lb::OriginMcRec,
                  hf_cand_lb::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandLbMcGen, "AOD", "HFCANDLBMCGEN", //!
                  hf_cand_lb::FlagMcMatchGen,
                  hf_cand_lb::OriginMcGen);

// specific B0 candidate properties
namespace hf_cand_b0
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // D index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level

// mapping of decay types
enum DecayType { B0ToDPi };

// B0(B0bar) → D∓ π±
template <typename T>
auto ctB0(const T& candidate)
{
  return candidate.ct(RecoDecay::getMassPDG(pdg::Code::kB0));
}

template <typename T>
auto yB0(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kB0));
}

template <typename T>
auto eB0(const T& candidate)
{
  return candidate.e(RecoDecay::getMassPDG(pdg::Code::kB0));
}

template <typename T>
auto invMassB0ToDPi(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(pdg::Code::kDMinus), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto cosThetaStarB0(const T& candidate)
{
  return candidate.cosThetaStar(array{RecoDecay::getMassPDG(pdg::Code::kDMinus), RecoDecay::getMassPDG(kPiPlus)}, RecoDecay::getMassPDG(pdg::Code::kB0), 1);
}
} // namespace hf_cand_b0

// declare dedicated B0 decay candidate table
DECLARE_SOA_TABLE(HfCandB0Base, "AOD", "HFCANDB0BASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand_b0::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::ImpactParameterProngSqSum<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandB0Ext, HfCandB0Base, "HFCANDB0EXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandB0 = HfCandB0Ext;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandB0McRec, "AOD", "HFCANDB0MCREC",
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandB0McGen, "AOD", "HFCANDB0MCGEN",
                  hf_cand_b0::FlagMcMatchGen,
                  hf_cand_b0::OriginMcGen);

// specific Σc0,++ candidate properties
namespace hf_cand_sigmac
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngLc, prongLc, int, HfCand3Prong, "");                //! Index to a Lc prong
DECLARE_SOA_COLUMN(Charge, charge, int8_t);                                            //! // Σc charge(either 0 or ++)
DECLARE_SOA_COLUMN(StatusSpreadLcMinvPKPiFromPDG, statusSpreadLcMinvPKPiFromPDG, int); //! // Λc Minv(pKpi) spread from PDG Λc mass
DECLARE_SOA_COLUMN(StatusSpreadLcMinvPiKPFromPDG, statusSpreadLcMinvPiKPFromPDG, int); //! // Λc Minv(piKp) spread from PDG Λc mass
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0");                //! Λc index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level

// mapping of decay types
enum DecayType { Sc0ToPKPiPi = 0,
                 ScplusplusToPKPiPi };

/// Σc0,++ → Λc+(→pK-π+) π-,+
/// @brief Sc inv. mass using reco mass for Lc in pKpi and PDG mass for pion
template <typename T, typename U>
auto invMassScRecoLcToPKPi(const T& candidateSc, const U& candidateLc)
{
  return candidateSc.m(array{static_cast<double>(hf_cand_3prong::invMassLcToPKPi(candidateLc)), RecoDecay::getMassPDG(kPiPlus)});
}

/// @brief Sc inv. mass using reco mass for Lc in piKp and PDG mass for pion
template <typename T, typename U>
auto invMassScRecoLcToPiKP(const T& candidateSc, const U& candidateLc)
{
  return candidateSc.m(array{static_cast<double>(hf_cand_3prong::invMassLcToPiKP(candidateLc)), RecoDecay::getMassPDG(kPiPlus)});
}

template <typename T>
auto ySc0(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kSigmaC0));
}

template <typename T>
auto yScPlusPlus(const T& candidate)
{
  return candidate.y(RecoDecay::getMassPDG(pdg::Code::kSigmaCPlusPlus));
}

} // namespace hf_cand_sigmac

// declare dedicated Σc0,++ decay candidate table
// NB: no topology for Σc0,++ (strong decay)
DECLARE_SOA_TABLE(HfCandScBase, "AOD", "HFCANDSCBASE",
                  o2::soa::Index<>,
                  // general columns
                  hf_cand::CollisionId,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  // hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  // hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  // hf_track_index::ProngLcId, hf_track_index::Prong1Id,
                  hf_cand_sigmac::ProngLcId, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* Σc0,++ specific columns */
                  hf_cand_sigmac::Charge,
                  hf_cand_sigmac::StatusSpreadLcMinvPKPiFromPDG, hf_cand_sigmac::StatusSpreadLcMinvPiKPFromPDG,
                  /* prong 0 */
                  // hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  /* prong 1 */
                  // hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_2prong::M2<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Pt2<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::P2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandScExt, HfCandScBase, "HFCANDSCEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);
using HfCandSc = HfCandScExt;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandScMcRec, "AOD", "HFCANDSCMCREC", //!
                  hf_cand_sigmac::FlagMcMatchRec,
                  hf_cand_sigmac::OriginMcRec);

// table with results of generation level MC matching
DECLARE_SOA_TABLE(HfCandScMcGen, "AOD", "HFCANDSCMCGEN", //!
                  hf_cand_sigmac::FlagMcMatchGen,
                  hf_cand_sigmac::OriginMcGen);

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_
