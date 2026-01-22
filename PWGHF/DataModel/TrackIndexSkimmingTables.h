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

/// \file TrackIndexSkimmingTables.h
/// \brief Definitions of tables produced by the candidate skimming
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef PWGHF_DATAMODEL_TRACKINDEXSKIMMINGTABLES_H_
#define PWGHF_DATAMODEL_TRACKINDEXSKIMMINGTABLES_H_

#include "PWGHF/Utils/utilsEvSelHf.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cstdint>
#include <vector>

namespace o2::aod
{
namespace hf_sel_collision
{

// ================
// Collision selection table
// ================

DECLARE_SOA_COLUMN(WhyRejectColl, whyRejectColl, o2::hf_evsel::HfCollisionRejectionMask); //!
} // namespace hf_sel_collision

DECLARE_SOA_TABLE(HfSelCollision, "AOD", "HFSELCOLLISION", //!
                  hf_sel_collision::WhyRejectColl);

// ================
// Track selection tables
// ================

namespace hf_sel_track
{
DECLARE_SOA_COLUMN(IsSelProng, isSelProng, uint32_t);           //!
DECLARE_SOA_COLUMN(IsIdentifiedPid, isIdentifiedPid, uint32_t); //!
DECLARE_SOA_COLUMN(IsPositive, isPositive, bool);               //!
} // namespace hf_sel_track

DECLARE_SOA_TABLE(HfSelTrack, "AOD", "HFSELTRACK", //!
                  hf_sel_track::IsSelProng,
                  hf_sel_track::IsIdentifiedPid,
                  hf_sel_track::IsPositive);

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

// ================
// Track index skim tables
// ================

namespace hf_track_index
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, Tracks, "_0"); //! Index to first prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, Tracks, "_1"); //! Index to second prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, Tracks, "_2"); //! Index to third prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong3, prong3, int, Tracks, "_3"); //! Index to fourth prong
DECLARE_SOA_INDEX_COLUMN_FULL(Prong4, prong4, int, Tracks, "_4"); //! Index to fifth prong
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                 //! Index to V0 prong
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);                       //! Index to cascade prong
DECLARE_SOA_COLUMN(HFflag, hfflag, uint8_t);                      //! Bitmap to store selection results, o2-linter: disable=name/o2-column (written to disk)

DECLARE_SOA_COLUMN(FlagD0ToKPi, flagD0ToKPi, uint8_t);           //!
DECLARE_SOA_COLUMN(FlagJpsiToEE, flagJpsiToEE, uint8_t);         //!
DECLARE_SOA_COLUMN(FlagJpsiToMuMu, flagJpsiToMuMu, uint8_t);     //!
DECLARE_SOA_COLUMN(FlagDplusToPiKPi, flagDplusToPiKPi, uint8_t); //!
DECLARE_SOA_COLUMN(FlagLcToPKPi, flagLcToPKPi, uint8_t);         //!
DECLARE_SOA_COLUMN(FlagDsToKKPi, flagDsToKKPi, uint8_t);         //!
DECLARE_SOA_COLUMN(FlagXicToPKPi, flagXicToPKPi, uint8_t);       //!
DECLARE_SOA_COLUMN(FlagDstarToD0Pi, flagDstarToD0Pi, uint8_t);   //!

DECLARE_SOA_COLUMN(MlProbSkimD0ToKPi, mlProbSkimD0ToKPi, std::vector<float>);           //! ML probabilities (background, prompt, non-prompt) for D0->Kpi
DECLARE_SOA_COLUMN(MlProbSkimDplusToPiKPi, mlProbSkimDplusToPiKPi, std::vector<float>); //! ML probabilities (background, prompt, non-prompt) for D+->Kpipi
DECLARE_SOA_COLUMN(MlProbSkimDsToKKPi, mlProbSkimDsToKKPi, std::vector<float>);         //! ML probabilities (background, prompt, non-prompt) for Ds->KKpi
DECLARE_SOA_COLUMN(MlProbSkimLcToPKPi, mlProbSkimLcToPKPi, std::vector<float>);         //! ML probabilities (background, prompt, non-prompt) for Lc->pKpi
DECLARE_SOA_COLUMN(MlProbSkimXicToPKPi, mlProbSkimXicToPKPi, std::vector<float>);       //! ML probabilities (background, prompt, non-prompt) for Xic->pKpi
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
                  hf_track_index::CollisionId,
                  hf_track_index::CascadeId,
                  hf_track_index::Prong0Id,
                  hf_track_index::HFflag);

using HfCascLf2Prong = HfCascLf2Prongs::iterator;

DECLARE_SOA_TABLE(HfCascLf3Prongs, "AOD", "HFCASCLF3PRONG", //! Table for HF 3 prong candidates with a Cascade
                  o2::soa::Index<>,
                  hf_track_index::CollisionId,
                  hf_track_index::CascadeId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);

using HfCascLf3Prong = HfCascLf3Prongs::iterator;

namespace hf_track_index
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngD0, prongD0, int, Hf2Prongs, ""); //! Index to a D0 prong
} // namespace hf_track_index

DECLARE_SOA_TABLE(HfDstars_000, "AOD", "HFDSTAR", //! D* -> D0pi candidates (Run 2 converted format)
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::ProngD0Id);

DECLARE_SOA_TABLE_VERSIONED(HfDstars_001, "AOD", "HFDSTAR", 1, //! D* -> D0pi candidates (Run 3 format)
                            o2::soa::Index<>,
                            hf_track_index::CollisionId,
                            hf_track_index::Prong0Id,
                            hf_track_index::ProngD0Id);

using HfDstars = HfDstars_001;
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

DECLARE_SOA_TABLE(HfCutStatusDstar, "AOD", "HFCUTSTATUSDST", //!
                  hf_track_index::FlagDstarToD0Pi);

DECLARE_SOA_TABLE(Hf2ProngMlProbs, "AOD", "HF2PRONGMLPROB", //! Table for ML scores of HF 2 prong candidates
                  hf_track_index::MlProbSkimD0ToKPi);

DECLARE_SOA_TABLE(Hf3ProngMlProbs, "AOD", "HF3PRONGMLPROB", //! Table for ML scores of HF 3 prong candidates
                  hf_track_index::MlProbSkimDplusToPiKPi,
                  hf_track_index::MlProbSkimLcToPKPi,
                  hf_track_index::MlProbSkimDsToKKPi,
                  hf_track_index::MlProbSkimXicToPKPi);

// ================
// Primary-vertex refit tables
// ================

namespace hf_pv_refit
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
} // namespace hf_pv_refit

DECLARE_SOA_TABLE(HfPvRefit2Prong, "AOD", "HFPVREFIT2PRONG", //!
                  hf_pv_refit::PvRefitX,
                  hf_pv_refit::PvRefitY,
                  hf_pv_refit::PvRefitZ,
                  hf_pv_refit::PvRefitSigmaX2,
                  hf_pv_refit::PvRefitSigmaXY,
                  hf_pv_refit::PvRefitSigmaY2,
                  hf_pv_refit::PvRefitSigmaXZ,
                  hf_pv_refit::PvRefitSigmaYZ,
                  hf_pv_refit::PvRefitSigmaZ2);

DECLARE_SOA_TABLE(HfPvRefit3Prong, "AOD", "HFPVREFIT3PRONG", //!
                  hf_pv_refit::PvRefitX,
                  hf_pv_refit::PvRefitY,
                  hf_pv_refit::PvRefitZ,
                  hf_pv_refit::PvRefitSigmaX2,
                  hf_pv_refit::PvRefitSigmaXY,
                  hf_pv_refit::PvRefitSigmaY2,
                  hf_pv_refit::PvRefitSigmaXZ,
                  hf_pv_refit::PvRefitSigmaYZ,
                  hf_pv_refit::PvRefitSigmaZ2,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfPvRefitDstar, "AOD", "HFPVREFITDSTAR", //!
                  hf_pv_refit::PvRefitX,
                  hf_pv_refit::PvRefitY,
                  hf_pv_refit::PvRefitZ,
                  hf_pv_refit::PvRefitSigmaX2,
                  hf_pv_refit::PvRefitSigmaXY,
                  hf_pv_refit::PvRefitSigmaY2,
                  hf_pv_refit::PvRefitSigmaXZ,
                  hf_pv_refit::PvRefitSigmaYZ,
                  hf_pv_refit::PvRefitSigmaZ2,
                  o2::soa::Marker<2>);

// ================
// Decay types stored in HFflag
// ================

namespace hf_cand_2prong
{
enum DecayType {
  D0ToPiK = 0,
  JpsiToEE,
  JpsiToMuMu,
  N2ProngDecays
};
} // namespace hf_cand_2prong

namespace hf_cand_3prong
{
enum DecayType {
  DplusToPiKPi = 0,
  LcToPKPi,
  DsToKKPi,
  XicToPKPi,
  CdToDeKPi,
  CtToTrKPi,
  ChToHeKPi,
  N3ProngDecays
};
} // namespace hf_cand_3prong

namespace hf_cand_casc_lf
{
enum DecayType2Prong {
  XiczeroOmegaczeroToXiPi = 0,
  OmegaczeroToOmegaPi,
  OmegaczeroToOmegaK,
  N2ProngDecays
};

enum DecayType3Prong {
  XicplusToXiPiPi = 0,
  N3ProngDecays
};
} // namespace hf_cand_casc_lf

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_TRACKINDEXSKIMMINGTABLES_H_
