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

#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"

#include "ALICE3/DataModel/ECAL.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/Utils/utilsPid.h"

namespace o2::aod
{
// Table aliases

using TracksWCov = soa::Join<Tracks, TracksCov>;
using TracksWDca = soa::Join<Tracks, TracksDCA>;
using TracksWExtra = soa::Join<Tracks, TracksExtra>;
using TracksWCovDca = soa::Join<Tracks, TracksCov, TracksDCA>;
using TracksWCovExtra = soa::Join<Tracks, TracksCov, TracksExtra>;
using TracksWDcaExtra = soa::Join<Tracks, TracksDCA, TracksExtra>;
using TracksWCovDcaExtra = soa::Join<Tracks, TracksCov, TracksDCA, TracksExtra>;

using TracksWMc = soa::Join<Tracks, McTrackLabels>;

using TracksPidEl = soa::Join<aod::pidTPCFullEl, aod::pidTOFFullEl>;
using TracksPidMu = soa::Join<aod::pidTPCFullMu, aod::pidTOFFullMu>;
using TracksPidPi = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi>;
using TracksPidKa = soa::Join<aod::pidTPCFullKa, aod::pidTOFFullKa>;
using TracksPidPr = soa::Join<aod::pidTPCFullPr, aod::pidTOFFullPr>;

using TracksPidTinyEl = soa::Join<aod::pidTPCEl, aod::pidTOFEl>;
using TracksPidTinyMu = soa::Join<aod::pidTPCMu, aod::pidTOFMu>;
using TracksPidTinyPi = soa::Join<aod::pidTPCPi, aod::pidTOFPi>;
using TracksPidTinyKa = soa::Join<aod::pidTPCKa, aod::pidTOFKa>;
using TracksPidTinyPr = soa::Join<aod::pidTPCPr, aod::pidTOFPr>;

// namespace pid_tpc_tof_utils
// {
// /// Function to combine TPC and TOF NSigma (for ML purposes)
// /// \param tpcNSigma is the (binned) NSigma separation in TPC (if tiny = true)
// /// \param tofNSigma is the (binned) NSigma separation in TOF (if tiny = true)
// /// \return Node containing the combined NSigma of TPC and TOF
// template <bool tiny = false, typename T1>
// o2::framework::expressions::Node combineNSigma(const T1& tpcNSigma, const T1& tofNSigma)
// {
//   float defaultNSigmaTolerance = .1f;
//   float defaultNSigma = -999.f + defaultNSigmaTolerance; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

//   if constexpr (tiny) {
//     auto tpcBinWidth = 1.f * pidtpc_tiny::binning::bin_width;
//     auto tofBinWidth = 1.f * pidtof_tiny::binning::bin_width;

//     return o2::framework::expressions::ifnode((tpcNSigma * tpcBinWidth > defaultNSigma) && (tofNSigma * tofBinWidth > defaultNSigma), o2::framework::expressions::nsqrt(.5f * tpcNSigma * tpcNSigma * tpcBinWidth * tpcBinWidth + .5f * tofNSigma * tofNSigma * tofBinWidth * tofBinWidth), // TPC and TOF
//                                               o2::framework::expressions::ifnode(tpcNSigma * tpcBinWidth > defaultNSigma, o2::framework::expressions::nabs(tpcNSigma * tpcBinWidth),                                                                                                        // only TPC
//                                                                                  o2::framework::expressions::ifnode(tofNSigma * tofBinWidth > defaultNSigma, o2::framework::expressions::nabs(tofNSigma * tofBinWidth),                                                                     // only TOF
//                                                                                                                     1.f * tofNSigma * tofBinWidth)));                                                                                                                                       // no TPC nor TOF
//   }

//   return o2::framework::expressions::ifnode((tpcNSigma > defaultNSigma) && (tofNSigma > defaultNSigma), o2::framework::expressions::nsqrt(.5f * tpcNSigma * tpcNSigma + .5f * tofNSigma * tofNSigma), // TPC and TOF
//                                             o2::framework::expressions::ifnode(tpcNSigma > defaultNSigma, o2::framework::expressions::nabs(tpcNSigma),                                                // only TPC
//                                                                                o2::framework::expressions::ifnode(tofNSigma > defaultNSigma, o2::framework::expressions::nabs(tofNSigma),             // only TOF
//                                                                                                                   1.f * tofNSigma)));                                                                 // no TPC nor TOF
// }
// } // namespace pid_tpc_tof_utils

// namespace pid_tpc_tof_full
// {
// // Combined TPC and TOF NSigma
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, //! Combined NSigma separation with the TPC & TOF detectors for electron
//                               float, pid_tpc_tof_utils::combineNSigma(o2::aod::pidtpc::tpcNSigmaEl, o2::aod::pidtof::tofNSigmaEl));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, //! Combined NSigma separation with the TPC & TOF detectors for muon
//                               float, pid_tpc_tof_utils::combineNSigma(o2::aod::pidtpc::tpcNSigmaMu, o2::aod::pidtof::tofNSigmaMu));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, //! Combined NSigma separation with the TPC & TOF detectors for pion
//                               float, pid_tpc_tof_utils::combineNSigma(o2::aod::pidtpc::tpcNSigmaPi, o2::aod::pidtof::tofNSigmaPi));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, //! Combined NSigma separation with the TPC & TOF detectors for kaon
//                               float, pid_tpc_tof_utils::combineNSigma(o2::aod::pidtpc::tpcNSigmaKa, o2::aod::pidtof::tofNSigmaKa));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, //! Combined NSigma separation with the TPC & TOF detectors for proton
//                               float, pid_tpc_tof_utils::combineNSigma(o2::aod::pidtpc::tpcNSigmaPr, o2::aod::pidtof::tofNSigmaPr));
// } // namespace pid_tpc_tof_full

// namespace pid_tpc_tof_tiny
// {
// // Combined binned TPC and TOF NSigma
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, //! Combined binned NSigma separation with the TPC & TOF detectors for electron
//                               float, pid_tpc_tof_utils::combineNSigma<true>(o2::aod::pidtpc_tiny::tpcNSigmaStoreEl, o2::aod::pidtof_tiny::tofNSigmaStoreEl));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, //! Combined binned NSigma separation with the TPC & TOF detectors for muon
//                               float, pid_tpc_tof_utils::combineNSigma<true>(o2::aod::pidtpc_tiny::tpcNSigmaStoreMu, o2::aod::pidtof_tiny::tofNSigmaStoreMu));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, //! Combined binned NSigma separation with the TPC & TOF detectors for pion
//                               float, pid_tpc_tof_utils::combineNSigma<true>(o2::aod::pidtpc_tiny::tpcNSigmaStorePi, o2::aod::pidtof_tiny::tofNSigmaStorePi));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, //! Combined binned NSigma separation with the TPC & TOF detectors for kaon
//                               float, pid_tpc_tof_utils::combineNSigma<true>(o2::aod::pidtpc_tiny::tpcNSigmaStoreKa, o2::aod::pidtof_tiny::tofNSigmaStoreKa));
// DECLARE_SOA_EXPRESSION_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, //! Combined binned NSigma separation with the TPC & TOF detectors for proton
//                               float, pid_tpc_tof_utils::combineNSigma<true>(o2::aod::pidtpc_tiny::tpcNSigmaStorePr, o2::aod::pidtof_tiny::tofNSigmaStorePr));
// } // namespace pid_tpc_tof_tiny

// // Extension of per particle tables
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidElExt, TracksPidEl, "PIDELEXT", //! Table of the TPC & TOF Combined NSigma for electron
//                                 pid_tpc_tof_full::TpcTofNSigmaEl);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidMuExt, TracksPidMu, "PIDMUEXT", //! Table of the TPC & TOF Combined NSigma for muon
//                                 pid_tpc_tof_full::TpcTofNSigmaMu);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidPiExt, TracksPidPi, "PIDPIEXT", //! Table of the TPC & TOF Combined NSigma for pion
//                                 pid_tpc_tof_full::TpcTofNSigmaPi);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidKaExt, TracksPidKa, "PIDKAEXT", //! Table of the TPC & TOF Combined NSigma for kaon
//                                 pid_tpc_tof_full::TpcTofNSigmaKa);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidPrExt, TracksPidPr, "PIDPREXT", //! Table of the TPC & TOF Combined NSigma for proton
//                                 pid_tpc_tof_full::TpcTofNSigmaPr);

// // Extension of tiny size tables
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidTinyElExt, TracksPidTinyEl, "PIDTINYELEXT", //! Table of the TPC & TOF combined binned NSigma for electron
//                                 pid_tpc_tof_tiny::TpcTofNSigmaEl);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidTinyMuExt, TracksPidTinyMu, "PIDTINYMUEXT", //! Table of the TPC & TOF combined binned NSigma for muon
//                                 pid_tpc_tof_tiny::TpcTofNSigmaMu);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidTinyPiExt, TracksPidTinyPi, "PIDTINYPIEXT", //! Table of the TPC & TOF combined binned NSigma for pion
//                                 pid_tpc_tof_tiny::TpcTofNSigmaPi);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidTinyKaExt, TracksPidTinyKa, "PIDTINYKAEXT", //! Table of the TPC & TOF combined binned NSigma for kaon
//                                 pid_tpc_tof_tiny::TpcTofNSigmaKa);
// DECLARE_SOA_EXTENDED_TABLE_USER(TracksPidTinyPrExt, TracksPidTinyPr, "PIDTINYPREXT", //! Table of the TPC & TOF combined binned NSigma for proton
//                                 pid_tpc_tof_tiny::TpcTofNSigmaPr);

namespace pid_tpc_tof_static_full
{
// Combined TPC and TOF NSigma
DECLARE_SOA_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, float); //! Combined NSigma separation with the TPC & TOF detectors for electron
DECLARE_SOA_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, float); //! Combined NSigma separation with the TPC & TOF detectors for muon
DECLARE_SOA_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, float); //! Combined NSigma separation with the TPC & TOF detectors for pion
DECLARE_SOA_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, float); //! Combined NSigma separation with the TPC & TOF detectors for kaon
DECLARE_SOA_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, float); //! Combined NSigma separation with the TPC & TOF detectors for proton
} // namespace pid_tpc_tof_static_full

namespace pid_tpc_tof_static_tiny
{
// Combined TPC and TOF NSigma
DECLARE_SOA_COLUMN(TpcTofNSigmaEl, tpcTofNSigmaEl, float); //! Combined NSigma separation with the TPC & TOF detectors for electron
DECLARE_SOA_COLUMN(TpcTofNSigmaMu, tpcTofNSigmaMu, float); //! Combined NSigma separation with the TPC & TOF detectors for muon
DECLARE_SOA_COLUMN(TpcTofNSigmaPi, tpcTofNSigmaPi, float); //! Combined NSigma separation with the TPC & TOF detectors for pion
DECLARE_SOA_COLUMN(TpcTofNSigmaKa, tpcTofNSigmaKa, float); //! Combined NSigma separation with the TPC & TOF detectors for kaon
DECLARE_SOA_COLUMN(TpcTofNSigmaPr, tpcTofNSigmaPr, float); //! Combined NSigma separation with the TPC & TOF detectors for proton
} // namespace pid_tpc_tof_static_tiny

// Extension of per particle tables
DECLARE_SOA_TABLE(PidTpcTofFullEl, "AOD", "PIDTPCTOFFULLEL", //! Table of the TPC & TOF Combined NSigma for electron
                  pid_tpc_tof_static_full::TpcTofNSigmaEl);
DECLARE_SOA_TABLE(PidTpcTofFullMu, "AOD", "PIDTPCTOFFULLMU", //! Table of the TPC & TOF Combined NSigma for muon
                  pid_tpc_tof_static_full::TpcTofNSigmaMu);
DECLARE_SOA_TABLE(PidTpcTofFullPi, "AOD", "PIDTPCTOFFULLPI", //! Table of the TPC & TOF Combined NSigma for pion
                  pid_tpc_tof_static_full::TpcTofNSigmaPi);
DECLARE_SOA_TABLE(PidTpcTofFullKa, "AOD", "PIDTPCTOFFULLKA", //! Table of the TPC & TOF Combined NSigma for kaon
                  pid_tpc_tof_static_full::TpcTofNSigmaKa);
DECLARE_SOA_TABLE(PidTpcTofFullPr, "AOD", "PIDTPCTOFFULLPR", //! Table of the TPC & TOF Combined NSigma for proton
                  pid_tpc_tof_static_full::TpcTofNSigmaPr);

// Extension of per particle tables
DECLARE_SOA_TABLE(PidTpcTofTinyEl, "AOD", "PIDTPCTOFTINYEL", //! Table of the TPC & TOF Combined NSigma for electron
                  pid_tpc_tof_static_tiny::TpcTofNSigmaEl);
DECLARE_SOA_TABLE(PidTpcTofTinyMu, "AOD", "PIDTPCTOFTINYMU", //! Table of the TPC & TOF Combined NSigma for muon
                  pid_tpc_tof_static_tiny::TpcTofNSigmaMu);
DECLARE_SOA_TABLE(PidTpcTofTinyPi, "AOD", "PIDTPCTOFTINYPI", //! Table of the TPC & TOF Combined NSigma for pion
                  pid_tpc_tof_static_tiny::TpcTofNSigmaPi);
DECLARE_SOA_TABLE(PidTpcTofTinyKa, "AOD", "PIDTPCTOFTINYKA", //! Table of the TPC & TOF Combined NSigma for kaon
                  pid_tpc_tof_static_tiny::TpcTofNSigmaKa);
DECLARE_SOA_TABLE(PidTpcTofTinyPr, "AOD", "PIDTPCTOFTINYPR", //! Table of the TPC & TOF Combined NSigma for proton
                  pid_tpc_tof_static_tiny::TpcTofNSigmaPr);

namespace hf_sel_collision
{
DECLARE_SOA_COLUMN(WhyRejectColl, whyRejectColl, uint16_t); //!
} // namespace hf_sel_collision

DECLARE_SOA_TABLE(HfSelCollision, "AOD", "HFSELCOLLISION", //!
                  hf_sel_collision::WhyRejectColl);

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

DECLARE_SOA_COLUMN(FlagDstarToD0Pi, flagDstarToD0Pi, uint8_t); //!

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
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter0, errorImpactParameter0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ0, impactParameterZ0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ0, errorImpactParameterZ0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised0, impactParameterZNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong1, pt2Prong1, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng1, pVectorProng1, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter1, errorImpactParameter1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ1, impactParameterZ1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ1, errorImpactParameterZ1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised1, impactParameterZNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float); //!
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float); //!
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng2, ptProng2, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong2, pt2Prong2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng2, pVectorProng2, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameter2, impactParameter2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameter2, errorImpactParameter2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ2, impactParameterZ2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ2, errorImpactParameterZ2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised2, impactParameterZNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
/// prong PID nsigma
DECLARE_SOA_COLUMN(NProngsContributorsPV, nProngsContributorsPV, uint8_t);           //! number of prongs contributing to the primary-vertex reconstruction
DECLARE_SOA_COLUMN(BitmapProngsContributorsPV, bitmapProngsContributorsPV, uint8_t); //! bitmap with booleans indicating prongs contributing to the primary-vertex reconstruction
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);                                   //! TPC nSigma for pion hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);                                   //! TPC nSigma for pion hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);                                   //! TPC nSigma for pion hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);                                   //! TPC nSigma for kaon hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);                                   //! TPC nSigma for kaon hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);                                   //! TPC nSigma for kaon hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);                                   //! TPC nSigma for proton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTpcPr1, nSigTpcPr1, float);                                   //! TPC nSigma for proton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);                                   //! TPC nSigma for proton hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);                                   //! TOF nSigma for pion hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);                                   //! TOF nSigma for pion hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);                                   //! TOF nSigma for pion hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);                                   //! TOF nSigma for kaon hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);                                   //! TOF nSigma for kaon hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);                                   //! TOF nSigma for kaon hypothesis - prong 2
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);                                   //! TOF nSigma for proton hypothesis - prong 0
DECLARE_SOA_COLUMN(NSigTofPr1, nSigTofPr1, float);                                   //! TOF nSigma for proton hypothesis - prong 1
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);                                   //! TOF nSigma for proton hypothesis - prong 2
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi0, tpcTofNSigmaPi0,                         //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 0
                           [](float tpcNSigmaPi0, float tofNSigmaPi0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi0, tofNSigmaPi0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi1, tpcTofNSigmaPi1, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 1
                           [](float tpcNSigmaPi1, float tofNSigmaPi1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi1, tofNSigmaPi1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPi2, tpcTofNSigmaPi2, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 2
                           [](float tpcNSigmaPi2, float tofNSigmaPi2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPi2, tofNSigmaPi2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa0, tpcTofNSigmaKa0, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 0
                           [](float tpcNSigmaKa0, float tofNSigmaKa0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa0, tofNSigmaKa0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa1, tpcTofNSigmaKa1, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 1
                           [](float tpcNSigmaKa1, float tofNSigmaKa1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa1, tofNSigmaKa1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaKa2, tpcTofNSigmaKa2, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 2
                           [](float tpcNSigmaKa2, float tofNSigmaKa2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaKa2, tofNSigmaKa2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr0, tpcTofNSigmaPr0, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 0
                           [](float tpcNSigmaPr0, float tofNSigmaPr0) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr0, tofNSigmaPr0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr1, tpcTofNSigmaPr1, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 1
                           [](float tpcNSigmaPr1, float tofNSigmaPr1) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr1, tofNSigmaPr1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaPr2, tpcTofNSigmaPr2, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 2
                           [](float tpcNSigmaPr2, float tofNSigmaPr2) -> float { return pid_tpc_tof_utils::combineNSigma<false /*tiny*/>(tpcNSigmaPr2, tofNSigmaPr2); });
// tiny (binned) option
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi0, tpcTofNSigmaTinyPi0, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 0
                           [](float tpcNSigmaPi0, float tofNSigmaPi0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi0, tofNSigmaPi0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi1, tpcTofNSigmaTinyPi1, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 1
                           [](float tpcNSigmaPi1, float tofNSigmaPi1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi1, tofNSigmaPi1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPi2, tpcTofNSigmaTinyPi2, //! Combined NSigma separation with the TPC & TOF detectors for pion - prong 2
                           [](float tpcNSigmaPi2, float tofNSigmaPi2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPi2, tofNSigmaPi2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa0, tpcTofNSigmaTinyKa0, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 0
                           [](float tpcNSigmaKa0, float tofNSigmaKa0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa0, tofNSigmaKa0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa1, tpcTofNSigmaTinyKa1, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 1
                           [](float tpcNSigmaKa1, float tofNSigmaKa1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa1, tofNSigmaKa1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyKa2, tpcTofNSigmaTinyKa2, //! Combined NSigma separation with the TPC & TOF detectors for kaon - prong 2
                           [](float tpcNSigmaKa2, float tofNSigmaKa2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaKa2, tofNSigmaKa2); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr0, tpcTofNSigmaTinyPr0, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 0
                           [](float tpcNSigmaPr0, float tofNSigmaPr0) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr0, tofNSigmaPr0); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr1, tpcTofNSigmaTinyPr1, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 1
                           [](float tpcNSigmaPr1, float tofNSigmaPr1) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr1, tofNSigmaPr1); });
DECLARE_SOA_DYNAMIC_COLUMN(TpcTofNSigmaTinyPr2, tpcTofNSigmaTinyPr2, //! Combined NSigma separation with the TPC & TOF detectors for proton - prong 2
                           [](float tpcNSigmaPr2, float tofNSigmaPr2) -> float { return pid_tpc_tof_utils::combineNSigma<true /*tiny*/>(tpcNSigmaPr2, tofNSigmaPr2); });

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
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //!
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::y(std::array{px, py, pz}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E2, e2, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e2(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);     //!
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa,                               //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CPAXY, cpaXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(Ct, ct, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_COLUMN(KfTopolChi2OverNdf, kfTopolChi2OverNdf, float); //! chi2overndf of the KFParticle topological constraint
// B-hadron mother information
DECLARE_SOA_COLUMN(PtBhadMotherPart, ptBhadMotherPart, float); //! pt of the first B-hadron mother particle (only in case of non-prompt)
DECLARE_SOA_COLUMN(PdgBhadMotherPart, pdgBhadMotherPart, int); //! pdg of the first B-hadron mother particle (only in case of non-prompt)
DECLARE_SOA_COLUMN(IdxBhadMotherPart, idxBhadMotherPart, int); //! index of the first B-hadron mother particle (only in case of non-prompt)
// Kink topology and material interaction mc flags
DECLARE_SOA_COLUMN(NTracksDecayed, nTracksDecayed, int8_t);                       //! number of tracks matched with kinked decay topology
DECLARE_SOA_COLUMN(NInteractionsWithMaterial, nInteractionsWithMaterial, int8_t); //! number of tracks matched after interaction with material

// method of secondary-vertex reconstruction
enum VertexerType { DCAFitter = 0,
                    KfParticle };
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
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStar, cosThetaStar, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m, double mTot, int iProng) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m, mTot, iProng); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1}, std::array{errIp0, errIp1}, std::array{std::array{px0, py0}, std::array{px1, py1}}); });
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level
// KF related properties
DECLARE_SOA_COLUMN(KfGeoMassD0, kfGeoMassD0, float);       //! mass of the D0 candidate from the KFParticle geometric fit
DECLARE_SOA_COLUMN(KfGeoMassD0bar, kfGeoMassD0bar, float); //! mass of the D0bar candidate from the KFParticle geometric fit

// mapping of decay types
enum DecayType { D0ToPiK = 0,
                 JpsiToEE,
                 JpsiToMuMu,
                 N2ProngDecays }; // always keep N2ProngDecays at the end

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
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_cand::NProngsContributorsPV, hf_cand::BitmapProngsContributorsPV,
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

DECLARE_SOA_TABLE(HfProng0PidPi, "AOD", "HFP0PIDPI", //!
                  hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0,
                  hf_cand::TpcTofNSigmaPi0<hf_cand::NSigTpcPi0, hf_cand::NSigTofPi0>);
DECLARE_SOA_TABLE(HfProng1PidPi, "AOD", "HFP1PIDPI", //!
                  hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1,
                  hf_cand::TpcTofNSigmaPi1<hf_cand::NSigTpcPi1, hf_cand::NSigTofPi1>);
DECLARE_SOA_TABLE(HfProng2PidPi, "AOD", "HFP2PIDPI", //!
                  hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2,
                  hf_cand::TpcTofNSigmaPi2<hf_cand::NSigTpcPi2, hf_cand::NSigTofPi2>);
DECLARE_SOA_TABLE(HfProng0PidKa, "AOD", "HFP0PIDKA", //!
                  hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0,
                  hf_cand::TpcTofNSigmaKa0<hf_cand::NSigTpcKa0, hf_cand::NSigTofKa0>);
DECLARE_SOA_TABLE(HfProng1PidKa, "AOD", "HFP1PIDKA", //!
                  hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1,
                  hf_cand::TpcTofNSigmaKa1<hf_cand::NSigTpcKa1, hf_cand::NSigTofKa1>);
DECLARE_SOA_TABLE(HfProng2PidKa, "AOD", "HFP2PIDKA", //!
                  hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2,
                  hf_cand::TpcTofNSigmaKa2<hf_cand::NSigTpcKa2, hf_cand::NSigTofKa2>);
DECLARE_SOA_TABLE(HfProng0PidPr, "AOD", "HFP0PIDPR", //!
                  hf_cand::NSigTpcPr0, hf_cand::NSigTofPr0,
                  hf_cand::TpcTofNSigmaPr0<hf_cand::NSigTpcPr0, hf_cand::NSigTofPr0>);
DECLARE_SOA_TABLE(HfProng1PidPr, "AOD", "HFP1PIDPR", //!
                  hf_cand::NSigTpcPr1, hf_cand::NSigTofPr1,
                  hf_cand::TpcTofNSigmaPr1<hf_cand::NSigTpcPr1, hf_cand::NSigTofPr1>);
DECLARE_SOA_TABLE(HfProng2PidPr, "AOD", "HFP2PIDPR", //!
                  hf_cand::NSigTpcPr2, hf_cand::NSigTofPr2,
                  hf_cand::TpcTofNSigmaPr2<hf_cand::NSigTpcPr2, hf_cand::NSigTofPr2>);

using HfCand2Prong = HfCand2ProngExt;
using HfCand2ProngWPid = soa::Join<HfCand2Prong, HfProng0PidPi, HfProng0PidKa, HfProng1PidPi, HfProng1PidKa>;

DECLARE_SOA_TABLE(HfCand2ProngKF, "AOD", "HFCAND2PKF",
                  hf_cand::KfTopolChi2OverNdf,
                  hf_cand_2prong::KfGeoMassD0, hf_cand_2prong::KfGeoMassD0bar);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcRec, "AOD", "HFCAND2PMCREC", //!
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand2ProngMcGen, "AOD", "HFCAND2PMCGEN", //!
                  hf_cand_2prong::FlagMcMatchGen,
                  hf_cand_2prong::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

// cascade decay candidate table

namespace hf_cand_casc
{
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! px of candidate
                              float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! py of candidate
                              float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! pz of candidate
                              float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
// DECLARE_SOA_DYNAMIC_COLUMN(M, m, [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Pos, ptV0Pos, //! pt of the positive V0 daughter
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0Neg, ptV0Neg, //! pt of the negative V0 daughter
                           [](float px, float py) { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(CtV0, ctV0, //! c*t of the V0
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level
DECLARE_SOA_COLUMN(V0X, v0x, float);                        //! X position of V0 decay
DECLARE_SOA_COLUMN(V0Y, v0y, float);                        //! Y position of V0 decay
DECLARE_SOA_COLUMN(V0Z, v0z, float);                        //! Z position of V0 decay
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
                  v0data::V0CosPA,
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
                  v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  hf_cand_casc::CtV0<hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_casc::V0X, hf_cand_casc::V0Y, hf_cand_casc::V0Z, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);
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
                  hf_cand_casc::FlagMcMatchRec,
                  hf_cand_casc::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandCascadeMcGen, "AOD", "HFCANDCASCMCGEN", //!
                  hf_cand_casc::FlagMcMatchGen,
                  hf_cand_casc::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

// specific BPlus candidate properties
namespace hf_cand_bplus
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand2Prong, "_0"); // D0 index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                 // debug flag for mis-association reconstruction level

enum DecayType { BplusToD0Pi = 0 };

enum DecayTypeMc : uint8_t { BplusToD0PiToKPiPi = 0,
                             PartlyRecoDecay,
                             OtherDecay,
                             NDecayTypeMc };
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

DECLARE_SOA_TABLE(HfCandBplusProngs, "AOD", "HFCANDBPPRONGS",
                  hf_cand_bplus::Prong0Id, hf_track_index::Prong1Id);

using HfCandBplus = soa::Join<HfCandBplusExt, HfCandBplusProngs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcRec, "AOD", "HFCANDBPMCREC",
                  hf_cand_bplus::FlagMcMatchRec,
                  hf_cand_bplus::OriginMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandBplusMcGen, "AOD", "HFCANDBPMCGEN",
                  hf_cand_bplus::FlagMcMatchGen,
                  hf_cand_bplus::OriginMcGen);

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
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const std::array<double, 3>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(M2, m2, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, float px2, float py2, float pz2, const std::array<double, 3>& m) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}, std::array{px2, py2, pz2}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSum, impactParameterProngSqSum, //!
                           [](float impParProng0, float impParProng1, float impParProng2) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1, impParProng2); });
DECLARE_SOA_DYNAMIC_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float ip2, float errIp2, float px0, float py0, float px1, float py1, float px2, float py2) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1, ip2}, std::array{errIp0, errIp1, errIp2}, std::array{std::array{px0, py0}, std::array{px1, py1}, std::array{px2, py2}}); });
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

// Ds± → K± K∓ π± or D± → K± K∓ π±

enum DecayChannelDToKKPi {
  DsToPhiPi = 1,
  DsToK0starK,
  DplusToPhiPi,  // used to describe D+ in MC production for Ds analysis
  DplusToK0starK // used to describe D+ in MC production for Ds analysis
};

// KF related properties
DECLARE_SOA_COLUMN(KfXPVError, kfXPVError, float);                       //! error of X coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfYPVError, kfYPVError, float);                       //! error of Y coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfZPVError, kfZPVError, float);                       //! error of Z coordinate of the event's primary vertex
DECLARE_SOA_COLUMN(KfXError, kfXError, float);                           //! error of candidate's decay point X coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfYError, kfYError, float);                           //! error of candidate's decay point Y coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfZError, kfZError, float);                           //! error of candidate's decay point Z coordinate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPKPi, kfMassPKPi, float);                       //! mass of the PKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKP, kfMassPiKP, float);                       //! mass of the PiKP candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKPi, kfMassPiKPi, float);                     //! mass of the PiKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassKKPi, kfMassKKPi, float);                       //! mass of the KKPi candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiKK, kfMassPiKK, float);                       //! mass of the PiKK candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassKPi, kfMassKPi, float);                         //! mass of the KPi pair from the KFParticle fit
DECLARE_SOA_COLUMN(KfMassPiK, kfMassPiK, float);                         //! mass of the PiK pair from the KFParticle fit
DECLARE_SOA_COLUMN(KfPx, kfPx, float);                                   //! Px of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfPy, kfPy, float);                                   //! Py of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfPz, kfPz, float);                                   //! Pz of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPx, kfErrorPx, float);                         //! Px error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPy, kfErrorPy, float);                         //! Py error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfErrorPz, kfErrorPz, float);                         //! Pz error of the candidate from the KFParticle fit
DECLARE_SOA_COLUMN(KfChi2PrimProng0, kfChi2PrimProng0, float);           //! chi2 primary of the first prong
DECLARE_SOA_COLUMN(KfChi2PrimProng1, kfChi2PrimProng1, float);           //! chi2 primary of the second prong
DECLARE_SOA_COLUMN(KfChi2PrimProng2, kfChi2PrimProng2, float);           //! chi2 primary of the third prong
DECLARE_SOA_COLUMN(KfDcaProng0Prong1, kfDcaProng0Prong1, float);         //! DCA between first and second prongs
DECLARE_SOA_COLUMN(KfDcaProng0Prong2, kfDcaProng0Prong2, float);         //! DCA between first and third prongs
DECLARE_SOA_COLUMN(KfDcaProng1Prong2, kfDcaProng1Prong2, float);         //! DCA between second and third prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng0Prong1, kfChi2GeoProng0Prong1, float); //! chi2 geo between first and second prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng0Prong2, kfChi2GeoProng0Prong2, float); //! chi2 geo between first and third prongs
DECLARE_SOA_COLUMN(KfChi2GeoProng1Prong2, kfChi2GeoProng1Prong2, float); //! chi2 geo between second and third prongs
DECLARE_SOA_COLUMN(KfChi2Geo, kfChi2Geo, float);                         //! chi2 geo of the full candidate
DECLARE_SOA_COLUMN(KfChi2Topo, kfChi2Topo, float);                       //! chi2 topo of the full candidate (chi2prim of candidate to PV)
DECLARE_SOA_COLUMN(KfDecayLength, kfDecayLength, float);                 //! decay length
DECLARE_SOA_COLUMN(KfDecayLengthError, kfDecayLengthError, float);       //! decay length error

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
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1, hf_cand::ImpactParameterZ2,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1, hf_cand::ErrorImpactParameterZ2,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id, hf_track_index::Prong2Id, hf_cand::NProngsContributorsPV, hf_cand::BitmapProngsContributorsPV,
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

DECLARE_SOA_TABLE(HfCand3ProngKF, "AOD", "HFCAND3PKF",
                  hf_cand_3prong::KfXError, hf_cand_3prong::KfYError, hf_cand_3prong::KfZError,
                  hf_cand_3prong::KfXPVError, hf_cand_3prong::KfYPVError, hf_cand_3prong::KfZPVError,
                  hf_cand_3prong::KfMassPKPi, hf_cand_3prong::KfMassPiKP, hf_cand_3prong::KfMassPiKPi, hf_cand_3prong::KfMassKKPi, hf_cand_3prong::KfMassPiKK, hf_cand_3prong::KfMassKPi, hf_cand_3prong::KfMassPiK,
                  hf_cand_3prong::KfPx, hf_cand_3prong::KfPy, hf_cand_3prong::KfPz,
                  hf_cand_3prong::KfErrorPx, hf_cand_3prong::KfErrorPy, hf_cand_3prong::KfErrorPz,
                  hf_cand_3prong::KfChi2PrimProng0, hf_cand_3prong::KfChi2PrimProng1, hf_cand_3prong::KfChi2PrimProng2,
                  hf_cand_3prong::KfDcaProng1Prong2, hf_cand_3prong::KfDcaProng0Prong2, hf_cand_3prong::KfDcaProng0Prong1,
                  hf_cand_3prong::KfChi2GeoProng1Prong2, hf_cand_3prong::KfChi2GeoProng0Prong2, hf_cand_3prong::KfChi2GeoProng0Prong1,
                  hf_cand_3prong::KfChi2Geo, hf_cand_3prong::KfDecayLength, hf_cand_3prong::KfDecayLengthError, hf_cand_3prong::KfChi2Topo);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcRec, "AOD", "HFCAND3PMCREC", //!
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::IsCandidateSwapped,
                  hf_cand_3prong::FlagMcDecayChanRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCand3ProngMcGen, "AOD", "HFCAND3PMCGEN", //!
                  hf_cand_3prong::FlagMcMatchGen,
                  hf_cand_3prong::OriginMcGen,
                  hf_cand_3prong::FlagMcDecayChanGen,
                  hf_cand::IdxBhadMotherPart);

namespace hf_cand_casc_lf
{
// mapping of decay types
enum DecayType2Prong { XiczeroOmegaczeroToXiPi = 0,
                       OmegaczeroToOmegaPi,
                       OmegaczeroToOmegaK,
                       N2ProngDecays }; // always keep N2ProngDecays at the end
// mapping of construct method
enum ConstructMethod { DcaFitter = 0,
                       KfParticle };
// mapping of decay types
enum DecayType3Prong { XicplusToXiPiPi = 0,
                       N3ProngDecays }; // always keep N3ProngDecays at the end
} // namespace hf_cand_casc_lf

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
namespace hf_cand_xic0_omegac0
{
// Data processing results:
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(XPv, xPv, float);
DECLARE_SOA_COLUMN(YPv, yPv, float);
DECLARE_SOA_COLUMN(ZPv, zPv, float);
DECLARE_SOA_COLUMN(XDecayVtxCharmBaryon, xDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(YDecayVtxCharmBaryon, yDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(ZDecayVtxCharmBaryon, zDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(XDecayVtxCascade, xDecayVtxCascade, float);
DECLARE_SOA_COLUMN(YDecayVtxCascade, yDecayVtxCascade, float);
DECLARE_SOA_COLUMN(ZDecayVtxCascade, zDecayVtxCascade, float);
DECLARE_SOA_COLUMN(XDecayVtxV0, xDecayVtxV0, float);
DECLARE_SOA_COLUMN(YDecayVtxV0, yDecayVtxV0, float);
DECLARE_SOA_COLUMN(ZDecayVtxV0, zDecayVtxV0, float);
DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t); // sign of pi <- xi
DECLARE_SOA_COLUMN(CovVtxCharmBaryon0, covVtxCharmBaryon0, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon1, covVtxCharmBaryon1, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon2, covVtxCharmBaryon2, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon3, covVtxCharmBaryon3, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon4, covVtxCharmBaryon4, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryon5, covVtxCharmBaryon5, float);
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PxBachFromCharmBaryon, pxBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyBachFromCharmBaryon, pyBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzBachFromCharmBaryon, pzBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PxLambda, pxLambda, float);
DECLARE_SOA_COLUMN(PyLambda, pyLambda, float);
DECLARE_SOA_COLUMN(PzLambda, pzLambda, float);
DECLARE_SOA_COLUMN(PxBachFromCasc, pxBachFromCasc, float);
DECLARE_SOA_COLUMN(PyBachFromCasc, pyBachFromCasc, float);
DECLARE_SOA_COLUMN(PzBachFromCasc, pzBachFromCasc, float);
DECLARE_SOA_COLUMN(PxPosV0Dau, pxPosV0Dau, float);
DECLARE_SOA_COLUMN(PyPosV0Dau, pyPosV0Dau, float);
DECLARE_SOA_COLUMN(PzPosV0Dau, pzPosV0Dau, float);
DECLARE_SOA_COLUMN(PxNegV0Dau, pxNegV0Dau, float);
DECLARE_SOA_COLUMN(PyNegV0Dau, pyNegV0Dau, float);
DECLARE_SOA_COLUMN(PzNegV0Dau, pzNegV0Dau, float);
DECLARE_SOA_COLUMN(ImpactParCascXY, impactParCascXY, float);
DECLARE_SOA_COLUMN(ImpactParBachFromCharmBaryonXY, impactParBachFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(ImpactParCascZ, impactParCascZ, float);
DECLARE_SOA_COLUMN(ImpactParBachFromCharmBaryonZ, impactParBachFromCharmBaryonZ, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParBachFromCharmBaryonXY, errImpactParBachFromCharmBaryonXY, float);
DECLARE_SOA_INDEX_COLUMN(V0, v0);
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);
DECLARE_SOA_INDEX_COLUMN_FULL(BachelorFromCharmBaryon, bachelorFromCharmBaryon, int, Tracks, "_bachelorfromcharmbaryon");
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(CosPACharmBaryon, cosPACharmBaryon, float);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, float);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, float);
DECLARE_SOA_COLUMN(CosPAXYCharmBaryon, cosPAXYCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, float);
DECLARE_SOA_COLUMN(CTauOmegac, ctauOmegac, float);
DECLARE_SOA_COLUMN(CTauCascade, ctauCascade, float);
DECLARE_SOA_COLUMN(CTauV0, ctauV0, float);
DECLARE_SOA_COLUMN(CTauXic, ctauXic, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaBachFromCasc, etaBachFromCasc, float);
DECLARE_SOA_COLUMN(EtaBachFromCharmBaryon, etaBachFromCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCharmBaryon, etaCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, float);
DECLARE_SOA_COLUMN(EtaV0, etaV0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float); // pos dau
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float); // neg dau
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau0, dcaZToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau1, dcaZToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaZToPvCascDau, dcaZToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(DecLenCharmBaryon, decLenCharmBaryon, float);
DECLARE_SOA_COLUMN(DecLenCascade, decLenCascade, float);
DECLARE_SOA_COLUMN(DecLenV0, decLenV0, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYCharmBaryon, errorDecayLengthXYCharmBaryon, float);

// KFParticle results
DECLARE_SOA_COLUMN(KfDcaXYPiFromOmegac, kfDcaXYPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfDcaXYPiFromXic, kfDcaXYPiFromXic, float);
DECLARE_SOA_COLUMN(KfDcaXYCascToPv, kfDcaXYCascToPv, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCasc, chi2GeoCasc, float);
DECLARE_SOA_COLUMN(Chi2GeoOmegac, chi2GeoOmegac, float);
DECLARE_SOA_COLUMN(Chi2GeoXic, chi2GeoXic, float);
DECLARE_SOA_COLUMN(Chi2MassV0, chi2MassV0, float);
DECLARE_SOA_COLUMN(Chi2MassCasc, chi2MassCasc, float);
DECLARE_SOA_COLUMN(V0ldl, v0ldl, float);
DECLARE_SOA_COLUMN(Cascldl, cascldl, float);
DECLARE_SOA_COLUMN(Omegacldl, omegacldl, float);
DECLARE_SOA_COLUMN(Xicldl, xicldl, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToPv, chi2TopoCascToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromOmegacToPv, chi2TopoPiFromOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromXicToPv, chi2TopoPiFromXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoOmegacToPv, chi2TopoOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoXicToPv, chi2TopoXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCasc, chi2TopoV0ToCasc, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToOmegac, chi2TopoCascToOmegac, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToXic, chi2TopoCascToXic, float);
DECLARE_SOA_COLUMN(DecayLenXYLambda, decayLenXYLambda, float);
DECLARE_SOA_COLUMN(DecayLenXYCasc, decayLenXYCasc, float);
DECLARE_SOA_COLUMN(DecayLenXYOmegac, decayLenXYOmegac, float);
DECLARE_SOA_COLUMN(DecayLenXYXic, decayLenXYXic, float);
DECLARE_SOA_COLUMN(CosPaV0ToCasc, cosPaV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaCascToOmegac, cosPaCascToOmegac, float);
DECLARE_SOA_COLUMN(CosPaCascToXic, cosPaCascToXic, float);
DECLARE_SOA_COLUMN(CosPaXYV0ToCasc, cosPaXYV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaXYCascToOmegac, cosPaXYCascToOmegac, float);
DECLARE_SOA_COLUMN(CosPaXYCascToXic, cosPaXYCascToXic, float);
DECLARE_SOA_COLUMN(KfMassV0, kfMassV0, float);
DECLARE_SOA_COLUMN(KfMassCasc, kfMassCasc, float);
DECLARE_SOA_COLUMN(KfMassOmegac, kfMassOmegac, float);
DECLARE_SOA_COLUMN(KfRapOmegac, kfRapOmegac, float);
DECLARE_SOA_COLUMN(KfRapXic, kfRapXic, float);
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfptPiFromXic, kfptPiFromXic, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(KfptXic, kfptXic, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromOmegac, cosThetaStarPiFromOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromXic, cosThetaStarPiFromXic, float);
DECLARE_SOA_COLUMN(V0Ndf, v0Ndf, float);
DECLARE_SOA_COLUMN(CascNdf, cascNdf, float);
DECLARE_SOA_COLUMN(OmegacNdf, omegacNdf, float);
DECLARE_SOA_COLUMN(XicNdf, xicNdf, float);
DECLARE_SOA_COLUMN(MassV0Ndf, massV0Ndf, float);
DECLARE_SOA_COLUMN(MassCascNdf, massCascNdf, float);
DECLARE_SOA_COLUMN(V0Chi2OverNdf, v0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(CascChi2OverNdf, cascChi2OverNdf, float);
DECLARE_SOA_COLUMN(OmegacChi2OverNdf, omegacChi2OverNdf, float);
DECLARE_SOA_COLUMN(XicChi2OverNdf, xicChi2OverNdf, float);
DECLARE_SOA_COLUMN(MassV0Chi2OverNdf, massV0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(MassCascChi2OverNdf, massCascChi2OverNdf, float);
DECLARE_SOA_COLUMN(CascRejectInvmass, cascRejectInvmass, float);

// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
DECLARE_SOA_COLUMN(DebugGenCharmBar, debugGenCharmBar, int8_t);
DECLARE_SOA_COLUMN(DebugGenCasc, debugGenCasc, int8_t);
DECLARE_SOA_COLUMN(DebugGenLambda, debugGenLambda, int8_t);
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(OriginGen, originGen, int8_t);
DECLARE_SOA_COLUMN(PtCharmBaryonGen, ptCharmBaryonGen, float);
DECLARE_SOA_COLUMN(RapidityCharmBaryonGen, rapidityCharmBaryonGen, float);

// dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(PtCharmBaryon, ptCharmBaryon,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtCasc, ptCasc,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtLambda, ptLambda,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtPiFromCharmBaryon, ptPiFromCharmBaryon,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtKaFromCasc, ptKaFromCasc,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });

// mapping of decay types
enum DecayType { XiczeroToXiPi = 0,
                 OmegaczeroToXiPi,
                 OmegaczeroToOmegaPi,
                 OmegaczeroToOmegaK };

} // end of namespace hf_cand_xic0_omegac0

// declare dedicated Omegac and Xic to Xi Pi candidate table
DECLARE_SOA_TABLE(HfCandToXiPi, "AOD", "HFCANDTOXIPI",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0, hf_cand_xic0_omegac0::CTauXic,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon);

DECLARE_SOA_TABLE(HfCandToOmegaPi, "AOD", "HFCANDTOOMEGAPI",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  // dynamic
                  hf_cand_xic0_omegac0::PtCharmBaryon<hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon>,
                  hf_cand_xic0_omegac0::PtCasc<hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc>,
                  hf_cand_xic0_omegac0::PtPiFromCharmBaryon<hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon>,
                  hf_cand_xic0_omegac0::PtLambda<hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda>,
                  hf_cand_xic0_omegac0::PtKaFromCasc<hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc>,

                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(HfCandToOmegaK, "AOD", "HFCANDTOOMEGAK",
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay, // charge pi<-cascade (neg -> omegac, pos -> antiomegac)
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon,
                  o2::soa::Marker<2>);

// table with results of KFParticle
DECLARE_SOA_TABLE(HfOmegacKf, "AOD", "HFOMEGACKF", //!
                  hf_cand_xic0_omegac0::KfDcaXYPiFromOmegac, hf_cand_xic0_omegac0::KfDcaXYCascToPv,
                  hf_cand_xic0_omegac0::Chi2GeoV0, hf_cand_xic0_omegac0::Chi2GeoCasc, hf_cand_xic0_omegac0::Chi2GeoOmegac,
                  hf_cand_xic0_omegac0::Chi2MassV0, hf_cand_xic0_omegac0::Chi2MassCasc,
                  hf_cand_xic0_omegac0::V0ldl, hf_cand_xic0_omegac0::Cascldl, hf_cand_xic0_omegac0::Omegacldl,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToPv, hf_cand_xic0_omegac0::Chi2TopoCascToPv, hf_cand_xic0_omegac0::Chi2TopoPiFromOmegacToPv, hf_cand_xic0_omegac0::Chi2TopoOmegacToPv,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToCasc, hf_cand_xic0_omegac0::Chi2TopoCascToOmegac,
                  hf_cand_xic0_omegac0::DecayLenXYLambda, hf_cand_xic0_omegac0::DecayLenXYCasc, hf_cand_xic0_omegac0::DecayLenXYOmegac,
                  hf_cand_xic0_omegac0::CosPaV0ToCasc, hf_cand_xic0_omegac0::CosPaCascToOmegac, hf_cand_xic0_omegac0::CosPaXYV0ToCasc, hf_cand_xic0_omegac0::CosPaXYCascToOmegac,
                  hf_cand_xic0_omegac0::KfRapOmegac,
                  hf_cand_xic0_omegac0::KfptPiFromOmegac, hf_cand_xic0_omegac0::KfptOmegac,
                  hf_cand_xic0_omegac0::CosThetaStarPiFromOmegac,
                  hf_cand_xic0_omegac0::V0Ndf, hf_cand_xic0_omegac0::CascNdf, hf_cand_xic0_omegac0::OmegacNdf,
                  hf_cand_xic0_omegac0::MassV0Ndf, hf_cand_xic0_omegac0::MassCascNdf,
                  hf_cand_xic0_omegac0::V0Chi2OverNdf, hf_cand_xic0_omegac0::CascChi2OverNdf, hf_cand_xic0_omegac0::OmegacChi2OverNdf,
                  hf_cand_xic0_omegac0::MassV0Chi2OverNdf, hf_cand_xic0_omegac0::MassCascChi2OverNdf, hf_cand_xic0_omegac0::CascRejectInvmass);

DECLARE_SOA_TABLE(HfCandToXiPiKf, "AOD", "HFCANDTOXIPIKF", //!
                  o2::soa::Index<>,
                  hf_cand_xic0_omegac0::CollisionId, hf_cand_xic0_omegac0::XPv, hf_cand_xic0_omegac0::YPv, hf_cand_xic0_omegac0::ZPv,
                  hf_cand_xic0_omegac0::XDecayVtxCharmBaryon, hf_cand_xic0_omegac0::YDecayVtxCharmBaryon, hf_cand_xic0_omegac0::ZDecayVtxCharmBaryon,
                  hf_cand_xic0_omegac0::XDecayVtxCascade, hf_cand_xic0_omegac0::YDecayVtxCascade, hf_cand_xic0_omegac0::ZDecayVtxCascade,
                  hf_cand_xic0_omegac0::XDecayVtxV0, hf_cand_xic0_omegac0::YDecayVtxV0, hf_cand_xic0_omegac0::ZDecayVtxV0,
                  hf_cand_xic0_omegac0::SignDecay,
                  hf_cand_xic0_omegac0::CovVtxCharmBaryon0, hf_cand_xic0_omegac0::CovVtxCharmBaryon1, hf_cand_xic0_omegac0::CovVtxCharmBaryon2, hf_cand_xic0_omegac0::CovVtxCharmBaryon3, hf_cand_xic0_omegac0::CovVtxCharmBaryon4, hf_cand_xic0_omegac0::CovVtxCharmBaryon5,
                  hf_cand_xic0_omegac0::PxCharmBaryon, hf_cand_xic0_omegac0::PyCharmBaryon, hf_cand_xic0_omegac0::PzCharmBaryon,
                  hf_cand_xic0_omegac0::PxCasc, hf_cand_xic0_omegac0::PyCasc, hf_cand_xic0_omegac0::PzCasc,
                  hf_cand_xic0_omegac0::PxBachFromCharmBaryon, hf_cand_xic0_omegac0::PyBachFromCharmBaryon, hf_cand_xic0_omegac0::PzBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::PxLambda, hf_cand_xic0_omegac0::PyLambda, hf_cand_xic0_omegac0::PzLambda,
                  hf_cand_xic0_omegac0::PxBachFromCasc, hf_cand_xic0_omegac0::PyBachFromCasc, hf_cand_xic0_omegac0::PzBachFromCasc,
                  hf_cand_xic0_omegac0::PxPosV0Dau, hf_cand_xic0_omegac0::PyPosV0Dau, hf_cand_xic0_omegac0::PzPosV0Dau,
                  hf_cand_xic0_omegac0::PxNegV0Dau, hf_cand_xic0_omegac0::PyNegV0Dau, hf_cand_xic0_omegac0::PzNegV0Dau,
                  hf_cand_xic0_omegac0::ImpactParCascXY, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonXY, hf_cand_xic0_omegac0::ImpactParCascZ, hf_cand_xic0_omegac0::ImpactParBachFromCharmBaryonZ,
                  hf_cand_xic0_omegac0::ErrImpactParCascXY, hf_cand_xic0_omegac0::ErrImpactParBachFromCharmBaryonXY,
                  hf_cand_xic0_omegac0::V0Id, v0data::PosTrackId, v0data::NegTrackId, hf_cand_xic0_omegac0::CascadeId, hf_cand_xic0_omegac0::BachelorFromCharmBaryonId, cascdata::BachelorId,
                  hf_cand_xic0_omegac0::InvMassLambda, hf_cand_xic0_omegac0::InvMassCascade, hf_cand_xic0_omegac0::InvMassCharmBaryon,
                  hf_cand_xic0_omegac0::CosPAV0, hf_cand_xic0_omegac0::CosPACharmBaryon, hf_cand_xic0_omegac0::CosPACasc, hf_cand_xic0_omegac0::CosPAXYV0, hf_cand_xic0_omegac0::CosPAXYCharmBaryon, hf_cand_xic0_omegac0::CosPAXYCasc,
                  hf_cand_xic0_omegac0::CTauOmegac, hf_cand_xic0_omegac0::CTauCascade, hf_cand_xic0_omegac0::CTauV0, hf_cand_xic0_omegac0::CTauXic,
                  hf_cand_xic0_omegac0::EtaV0PosDau, hf_cand_xic0_omegac0::EtaV0NegDau, hf_cand_xic0_omegac0::EtaBachFromCasc, hf_cand_xic0_omegac0::EtaBachFromCharmBaryon,
                  hf_cand_xic0_omegac0::EtaCharmBaryon, hf_cand_xic0_omegac0::EtaCascade, hf_cand_xic0_omegac0::EtaV0,
                  hf_cand_xic0_omegac0::DcaXYToPvV0Dau0, hf_cand_xic0_omegac0::DcaXYToPvV0Dau1, hf_cand_xic0_omegac0::DcaXYToPvCascDau,
                  hf_cand_xic0_omegac0::DcaZToPvV0Dau0, hf_cand_xic0_omegac0::DcaZToPvV0Dau1, hf_cand_xic0_omegac0::DcaZToPvCascDau,
                  hf_cand_xic0_omegac0::DcaCascDau, hf_cand_xic0_omegac0::DcaV0Dau, hf_cand_xic0_omegac0::DcaCharmBaryonDau,
                  hf_cand_xic0_omegac0::DecLenCharmBaryon, hf_cand_xic0_omegac0::DecLenCascade, hf_cand_xic0_omegac0::DecLenV0, hf_cand_xic0_omegac0::ErrorDecayLengthCharmBaryon, hf_cand_xic0_omegac0::ErrorDecayLengthXYCharmBaryon,
                  hf_cand_xic0_omegac0::KfDcaXYPiFromXic, hf_cand_xic0_omegac0::KfDcaXYCascToPv,
                  hf_cand_xic0_omegac0::Chi2GeoV0, hf_cand_xic0_omegac0::Chi2GeoCasc, hf_cand_xic0_omegac0::Chi2GeoXic,
                  hf_cand_xic0_omegac0::Chi2MassV0, hf_cand_xic0_omegac0::Chi2MassCasc,
                  hf_cand_xic0_omegac0::V0ldl, hf_cand_xic0_omegac0::Cascldl, hf_cand_xic0_omegac0::Xicldl,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToPv, hf_cand_xic0_omegac0::Chi2TopoCascToPv, hf_cand_xic0_omegac0::Chi2TopoPiFromXicToPv, hf_cand_xic0_omegac0::Chi2TopoXicToPv,
                  hf_cand_xic0_omegac0::Chi2TopoV0ToCasc, hf_cand_xic0_omegac0::Chi2TopoCascToXic,
                  hf_cand_xic0_omegac0::DecayLenXYLambda, hf_cand_xic0_omegac0::DecayLenXYCasc, hf_cand_xic0_omegac0::DecayLenXYXic,
                  hf_cand_xic0_omegac0::CosPaV0ToCasc, hf_cand_xic0_omegac0::CosPaCascToXic, hf_cand_xic0_omegac0::CosPaXYV0ToCasc, hf_cand_xic0_omegac0::CosPaXYCascToXic,
                  hf_cand_xic0_omegac0::KfRapXic,
                  hf_cand_xic0_omegac0::KfptPiFromXic, hf_cand_xic0_omegac0::KfptXic,
                  hf_cand_xic0_omegac0::CosThetaStarPiFromXic,
                  hf_cand_xic0_omegac0::V0Ndf, hf_cand_xic0_omegac0::CascNdf, hf_cand_xic0_omegac0::XicNdf,
                  hf_cand_xic0_omegac0::MassV0Ndf, hf_cand_xic0_omegac0::MassCascNdf,
                  hf_cand_xic0_omegac0::V0Chi2OverNdf, hf_cand_xic0_omegac0::CascChi2OverNdf, hf_cand_xic0_omegac0::XicChi2OverNdf,
                  hf_cand_xic0_omegac0::MassV0Chi2OverNdf, hf_cand_xic0_omegac0::MassCascChi2OverNdf);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfXicToXiPiMCRec, "AOD", "HFXICXIPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<1>);
DECLARE_SOA_TABLE(HfOmegacToXiPiMCRec, "AOD", "HFOMCXIPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<2>);
DECLARE_SOA_TABLE(HfToOmegaPiMCRec, "AOD", "HFTOOMEPIMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<3>);
DECLARE_SOA_TABLE(HfToOmegaKMCRec, "AOD", "HFTOOMEKMCREC", //!
                  hf_cand_xic0_omegac0::FlagMcMatchRec,
                  hf_cand_xic0_omegac0::DebugMcRec,
                  hf_cand_xic0_omegac0::OriginRec,
                  hf_cand_xic0_omegac0::CollisionMatched,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  o2::soa::Marker<4>);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfXicToXiPiMCGen, "AOD", "HFXICXIPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<1>);
DECLARE_SOA_TABLE(HfOmegacToXiPiMCGen, "AOD", "HFOMECXIPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<2>);
DECLARE_SOA_TABLE(HfToOmegaPiMCGen, "AOD", "HFTOOMEPIMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<3>);
DECLARE_SOA_TABLE(HfToOmegaKMCGen, "AOD", "HFTOOMEKMCGEN", //!
                  hf_cand_xic0_omegac0::FlagMcMatchGen, hf_cand_xic0_omegac0::DebugGenCharmBar, hf_cand_xic0_omegac0::DebugGenCasc, hf_cand_xic0_omegac0::DebugGenLambda,
                  hf_cand_xic0_omegac0::PtCharmBaryonGen, hf_cand_xic0_omegac0::RapidityCharmBaryonGen, hf_cand_xic0_omegac0::OriginGen, hf_cand::IdxBhadMotherPart, o2::soa::Marker<4>);

// specific Xic to Xi Pi Pi candidate properties
namespace hf_cand_xic_to_xi_pi_pi
{
DECLARE_SOA_INDEX_COLUMN_FULL(Pi0, pi0, int, Tracks, "_pi0");
DECLARE_SOA_INDEX_COLUMN_FULL(Pi1, pi1, int, Tracks, "_pi1");
DECLARE_SOA_COLUMN(XPvErr, xPvErr, float);
DECLARE_SOA_COLUMN(YPvErr, yPvErr, float);
DECLARE_SOA_COLUMN(ZPvErr, zPvErr, float);
DECLARE_SOA_COLUMN(XSvErr, xSvErr, float);
DECLARE_SOA_COLUMN(YSvErr, ySvErr, float);
DECLARE_SOA_COLUMN(ZSvErr, zSvErr, float);
DECLARE_SOA_COLUMN(CosPaXi, cosPaXi, float);
DECLARE_SOA_COLUMN(CosPaXYXi, cosPaXYXi, float);
DECLARE_SOA_COLUMN(CosPaLambda, cosPaLambda, float);
DECLARE_SOA_COLUMN(CosPaXYLambda, cosPaXYLambda, float);
DECLARE_SOA_COLUMN(CosPaLambdaToXi, cosPaLambdaToXi, float);
DECLARE_SOA_COLUMN(CosPaXYLambdaToXi, cosPaXYLambdaToXi, float);
DECLARE_SOA_COLUMN(InvMassXicPlus, invMassXicPlus, float);
DECLARE_SOA_COLUMN(InvMassXi, invMassXi, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(InvMassXiPi0, invMassXiPi0, float);
DECLARE_SOA_COLUMN(InvMassXiPi1, invMassXiPi1, float);
DECLARE_SOA_COLUMN(PBachelorPi, pBachelorPi, float);
DECLARE_SOA_COLUMN(PPiFromLambda, pPiFromLambda, float);
DECLARE_SOA_COLUMN(PPrFromLambda, pPrFromLambda, float);
DECLARE_SOA_COLUMN(DcaXiDaughters, dcaXiDaughters, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);
DECLARE_SOA_COLUMN(DcaBachelorToPV, dcaBachelorToPV, float);
DECLARE_SOA_COLUMN(DcaXYCascToPV, dcaXYCascToPV, float);
DECLARE_SOA_COLUMN(DcaZCascToPV, dcaZCascToPV, float);
// KF specific columns
DECLARE_SOA_COLUMN(DcaXYPi0Pi1, dcaXYPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaXYPi0Xi, dcaXYPi0Xi, float);
DECLARE_SOA_COLUMN(DcaXYPi1Xi, dcaXYPi1Xi, float);
DECLARE_SOA_COLUMN(DcaPi0Pi1, dcaPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaPi0Xi, dcaPi0Xi, float);
DECLARE_SOA_COLUMN(DcaPi1Xi, dcaPi1Xi, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPV, chi2TopoXicPlusToPV, float);
DECLARE_SOA_COLUMN(Chi2TopoXicPlusToPVBeforeConstraint, chi2TopoXicPlusToPVBeforeConstraint, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlus, chi2TopoXiToXicPlus, float);
DECLARE_SOA_COLUMN(Chi2TopoXiToXicPlusBeforeConstraint, chi2TopoXiToXicPlusBeforeConstraint, float);
// PID
DECLARE_SOA_COLUMN(NSigTpcPiFromXicPlus0, nSigTpcPiFromXicPlus0, float);
DECLARE_SOA_COLUMN(NSigTpcPiFromXicPlus1, nSigTpcPiFromXicPlus1, float);
DECLARE_SOA_COLUMN(NSigTpcBachelorPi, nSigTpcBachelorPi, float);
DECLARE_SOA_COLUMN(NSigTpcPiFromLambda, nSigTpcPiFromLambda, float);
DECLARE_SOA_COLUMN(NSigTpcPrFromLambda, nSigTpcPrFromLambda, float);
DECLARE_SOA_COLUMN(NSigTofPiFromXicPlus0, nSigTofPiFromXicPlus0, float);
DECLARE_SOA_COLUMN(NSigTofPiFromXicPlus1, nSigTofPiFromXicPlus1, float);
DECLARE_SOA_COLUMN(NSigTofBachelorPi, nSigTofBachelorPi, float);
DECLARE_SOA_COLUMN(NSigTofPiFromLambda, nSigTofPiFromLambda, float);
DECLARE_SOA_COLUMN(NSigTofPrFromLambda, nSigTofPrFromLambda, float);
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); // generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(DebugMcGen, debugMcGen, int8_t);
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(OriginGen, originGen, int8_t);
// Dynamic columns
DECLARE_SOA_DYNAMIC_COLUMN(PProng0, pProng0, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PProng1, pProng1, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PProng2, pProng2, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });

// mapping of decay types
enum DecayType { XicToXiPiPi = 0,      // Ξc± → Ξ∓ π± π±
                 XicToXiResPiToXiPiPi, // Ξc± → Ξ(1530) π± → Ξ∓ π± π±
                 NDecayType };
} // end of namespace hf_cand_xic_to_xi_pi_pi

// declare dedicated Xic to Xi Pi Pi candidate table
DECLARE_SOA_TABLE(HfCandXicBase, "AOD", "HFCANDXICBASE",
                  hf_cand::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_xic_to_xi_pi_pi::XPvErr, hf_cand_xic_to_xi_pi_pi::YPvErr, hf_cand_xic_to_xi_pi_pi::ZPvErr,
                  // 3-prong specific columns
                  cascdata::CascadeId, hf_cand_xic_to_xi_pi_pi::Pi0Id, hf_cand_xic_to_xi_pi_pi::Pi1Id,
                  cascdata::BachelorId, cascdata::PosTrackId, cascdata::NegTrackId,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand_xic_to_xi_pi_pi::XSvErr, hf_cand_xic_to_xi_pi_pi::YSvErr, hf_cand_xic_to_xi_pi_pi::ZSvErr,
                  hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA, hf_cand_xic_to_xi_pi_pi::InvMassXicPlus, hf_cand_xic_to_xi_pi_pi::Sign,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1, hf_cand::ErrorImpactParameter2,
                  // cascade specific columns
                  hf_cand_xic_to_xi_pi_pi::PBachelorPi, hf_cand_xic_to_xi_pi_pi::PPiFromLambda, hf_cand_xic_to_xi_pi_pi::PPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::CosPaXi, hf_cand_xic_to_xi_pi_pi::CosPaXYXi, hf_cand_xic_to_xi_pi_pi::CosPaLambda, hf_cand_xic_to_xi_pi_pi::CosPaXYLambda, hf_cand_xic_to_xi_pi_pi::CosPaLambdaToXi, hf_cand_xic_to_xi_pi_pi::CosPaXYLambdaToXi,
                  hf_cand_xic_to_xi_pi_pi::InvMassXi, hf_cand_xic_to_xi_pi_pi::InvMassLambda, hf_cand_xic_to_xi_pi_pi::InvMassXiPi0, hf_cand_xic_to_xi_pi_pi::InvMassXiPi1,
                  // DCA
                  hf_cand_xic_to_xi_pi_pi::DcaXiDaughters, hf_cand_xic_to_xi_pi_pi::DcaV0Daughters, hf_cand_xic_to_xi_pi_pi::DcaPosToPV, hf_cand_xic_to_xi_pi_pi::DcaNegToPV, hf_cand_xic_to_xi_pi_pi::DcaBachelorToPV, hf_cand_xic_to_xi_pi_pi::DcaXYCascToPV, hf_cand_xic_to_xi_pi_pi::DcaZCascToPV,
                  // PID
                  hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus0, hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromXicPlus1, hf_cand_xic_to_xi_pi_pi::NSigTpcBachelorPi, hf_cand_xic_to_xi_pi_pi::NSigTpcPiFromLambda, hf_cand_xic_to_xi_pi_pi::NSigTpcPrFromLambda,
                  hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus0, hf_cand_xic_to_xi_pi_pi::NSigTofPiFromXicPlus1, hf_cand_xic_to_xi_pi_pi::NSigTofBachelorPi, hf_cand_xic_to_xi_pi_pi::NSigTofPiFromLambda, hf_cand_xic_to_xi_pi_pi::NSigTofPrFromLambda,
                  /* dynamic columns */
                  hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,
                  hf_cand::DecayLengthXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ErrorDecayLength>,
                  hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY>,
                  hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  hf_cand::ImpactParameterNormalised2<hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2>,
                  /* dynamic columns that use daughter momentum components */
                  hf_cand_xic_to_xi_pi_pi::PProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand_xic_to_xi_pi_pi::PProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand_xic_to_xi_pi_pi::PProng2<hf_cand::PxProng2, hf_cand::PyProng2, hf_cand::PzProng2>,
                  hf_cand::PtProng2<hf_cand::PxProng2, hf_cand::PyProng2>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand::Pt<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::P<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::PVector<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::CPAXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand_3prong::MaxNormalisedDeltaIP<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY, hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::ImpactParameter2, hf_cand::ErrorImpactParameter2, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PxProng2, hf_cand::PyProng2>,
                  hf_cand::Eta<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>,
                  hf_cand::Phi<hf_cand_3prong::Px, hf_cand_3prong::Py>,
                  hf_cand::Y<hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandXicExt, HfCandXicBase, "HFCANDXICEXT",
                                hf_cand_3prong::Px, hf_cand_3prong::Py, hf_cand_3prong::Pz);

using HfCandXic = HfCandXicExt;

DECLARE_SOA_TABLE(HfCandXicKF, "AOD", "HFCANDXICKF",
                  cascdata::KFCascadeChi2, cascdata::KFV0Chi2,
                  hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPVBeforeConstraint, hf_cand_xic_to_xi_pi_pi::Chi2TopoXicPlusToPV, hf_cand_xic_to_xi_pi_pi::Chi2TopoXiToXicPlusBeforeConstraint, hf_cand_xic_to_xi_pi_pi::Chi2TopoXiToXicPlus,
                  hf_cand_xic_to_xi_pi_pi::DcaXYPi0Pi1, hf_cand_xic_to_xi_pi_pi::DcaXYPi0Xi, hf_cand_xic_to_xi_pi_pi::DcaXYPi1Xi,
                  hf_cand_xic_to_xi_pi_pi::DcaPi0Pi1, hf_cand_xic_to_xi_pi_pi::DcaPi0Xi, hf_cand_xic_to_xi_pi_pi::DcaPi1Xi,
                  cascdata::DCACascDaughters);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandXicMcRec, "AOD", "HFCANDXICMCREC", //!
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchRec,
                  hf_cand_xic_to_xi_pi_pi::DebugMcRec,
                  hf_cand_xic_to_xi_pi_pi::OriginRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandXicMcGen, "AOD", "HFCANDXICMCGEN", //!
                  hf_cand_xic_to_xi_pi_pi::FlagMcMatchGen,
                  hf_cand_xic_to_xi_pi_pi::DebugMcGen,
                  hf_cand_xic_to_xi_pi_pi::OriginGen);

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
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                 // debug flag for mis-association reconstruction level

// mapping of decay types
enum DecayType { B0ToDPi };

enum DecayTypeMc : uint8_t { B0ToDplusPiToPiKPiPi = 0,
                             B0ToDsPiToKKPiPi,
                             BsToDsPiToKKPiPi,
                             PartlyRecoDecay,
                             OtherDecay,
                             NDecayTypeMc };
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

DECLARE_SOA_TABLE(HfCandB0Prongs, "AOD", "HFCANDB0PRONGS",
                  hf_cand_b0::Prong0Id, hf_track_index::Prong1Id);

using HfCandB0 = soa::Join<HfCandB0Ext, HfCandB0Prongs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandB0McRec, "AOD", "HFCANDB0MCREC",
                  hf_cand_b0::FlagMcMatchRec,
                  hf_cand_b0::OriginMcRec,
                  hf_cand_b0::DebugMcRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandB0McGen, "AOD", "HFCANDB0MCGEN",
                  hf_cand_b0::FlagMcMatchGen,
                  hf_cand_b0::OriginMcGen);

// specific Bs candidate properties
namespace hf_cand_bs
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCand3Prong, "_0"); // Ds index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         // reconstruction level
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         // generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               // particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               // particle origin, generator level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                 // debug flag for mis-association reconstruction level

// mapping of decay types
enum DecayType { BsToDsPi };

enum DecayTypeMc : uint8_t { BsToDsPiToPhiPiPiToKKPiPi = 0, // Bs(bar) → Ds∓ π± → (Phi π∓) π± → (K- K+ π∓) π±
                             BsToDsPiToK0starKPiToKKPiPi,   // Bs(bar) → Ds∓ π± → (K0* K∓) π± → (K- K+ π∓) π±
                             B0ToDsPiToPhiPiPiToKKPiPi,     // B0(bar) → Ds± π∓ → (Phi π±) π∓ → (K- K+ π±) π∓
                             B0ToDsPiToK0starKPiToKKPiPi,   // B0(bar) → Ds± π∓ → (K0* K±) π∓ → (K- K+ π±) π∓
                             PartlyRecoDecay,               // 4 final state particles have another common b-hadron ancestor
                             OtherDecay,
                             NDecayTypeMc }; // counter of differentiated MC decay types

} // namespace hf_cand_bs

// declare dedicated Bs decay candidate table
DECLARE_SOA_TABLE(HfCandBsBase, "AOD", "HFCANDBSBASE",
                  // general columns
                  HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
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
                  hf_cand::E2<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  o2::soa::Marker<1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandBsExt, HfCandBsBase, "HFCANDBSEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

DECLARE_SOA_TABLE(HfCandBsProngs, "AOD", "HFCANDBSPRONGS",
                  hf_cand_bs::Prong0Id, hf_track_index::Prong1Id);

using HfCandBs = soa::Join<HfCandBsExt, HfCandBsProngs>;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandBsMcRec, "AOD", "HFCANDBSMCREC",
                  hf_cand_bs::FlagMcMatchRec);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandBsMcGen, "AOD", "HFCANDBSMCGEN",
                  hf_cand_bs::FlagMcMatchGen);

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
                  hf_cand_sigmac::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart);

// table with results of generation level MC matching
DECLARE_SOA_TABLE(HfCandScMcGen, "AOD", "HFCANDSCMCGEN", //!
                  hf_cand_sigmac::FlagMcMatchGen,
                  hf_cand_sigmac::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

// specific Σc0,++ candidate properties in cascade channel
namespace hf_cand_sigmac_to_cascade
{
DECLARE_SOA_INDEX_COLUMN_FULL(ProngLc, prongLc, int, HfCandCascade, "");               //! Index to a Lc prong
DECLARE_SOA_COLUMN(Charge, charge, int8_t);                                            //! // Σc charge(either 0 or ++)
DECLARE_SOA_COLUMN(ChargeLc, chargelc, int8_t);                                        //! // Λc charge(+)
DECLARE_SOA_COLUMN(ChargeSoftPi, chargesoftpi, int8_t);                                //! // pion charge(either - or +)
DECLARE_SOA_COLUMN(StatusSpreadLcMinvKs0PFromPDG, statusSpreadLcMinvKs0PFromPDG, int); //! // Λc Minv spread from PDG Λc mass
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, HfCandCascade, "_0");               //! Λc index
// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);       //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);       //! particle origin, generator level
} // namespace hf_cand_sigmac_to_cascade

// declare dedicated Σc0,++ decay candidate table
// NB: no topology for Σc0,++ (strong decay)
DECLARE_SOA_TABLE(HfCandScCasBase, "AOD", "HFCANDSCCASBASE",
                  o2::soa::Index<>,
                  // general columns
                  hf_cand::CollisionId,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand_sigmac_to_cascade::ProngLcId, hf_track_index::Prong1Id,
                  hf_cand_sigmac_to_cascade::ChargeLc,
                  hf_cand_sigmac_to_cascade::ChargeSoftPi,
                  // hf_track_index::HFflag,
                  /* Σc0,++ specific columns */
                  hf_cand_sigmac_to_cascade::Charge,
                  // hf_cand_sigmac_to_cascade::StatusSpreadLcMinvKs0PFromPDG,
                  /* prong 0 */
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  /* prong 1 */
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
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandScCasExt, HfCandScCasBase, "HFCANDSCCASEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);
using HfCandScCascades = HfCandScCasExt;
using HfCandScCascade = HfCandScCascades::iterator;

/// D*± → D0(bar) π±
namespace hf_cand_dstar
{
DECLARE_SOA_EXPRESSION_COLUMN(PxD0, pxD0, float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(PyD0, pyD0, float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(PzD0, pzD0, float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProductD0, impactParameterProductD0,
                           [](float dca1, float dca2) -> float { return dca1 * dca2; });
// Dynamic Columns for D0 candidate using PDG masses of daughters
DECLARE_SOA_DYNAMIC_COLUMN(InvMassD0, invMassD0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMass2D0, invMass2D0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarD0, cosThetaStarD0,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus}, constants::physics::MassD0, 1); });
// Dynamic Columns for D0Bar candidate using PDG masses of daughters
DECLARE_SOA_DYNAMIC_COLUMN(InvMassD0Bar, invMassD0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(InvMass2D0Bar, invMass2D0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::m2(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}); });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarD0Bar, cosThetaStarD0Bar,
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1) -> float { return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, std::array{constants::physics::MassKPlus, constants::physics::MassPiPlus}, constants::physics::MassD0, 0); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterProngSqSumD0, impactParameterProngSqSumD0,
                           [](float impParProng0, float impParProng1) -> float { return RecoDecay::sumOfSquares(impParProng0, impParProng1); });
DECLARE_SOA_DYNAMIC_COLUMN(DeltaIPNormalisedMaxD0, deltaIPNormalisedMaxD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float errDlxy, float pxM, float pyM, float ip0, float errIp0, float ip1, float errIp1, float px0, float py0, float px1, float py1) -> float { return RecoDecay::maxNormalisedDeltaIP(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, errDlxy, std::array{pxM, pyM}, std::array{ip0, ip1}, std::array{errIp0, errIp1}, std::array{std::array{px0, py0}, std::array{px1, py1}}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtD0, ptD0,
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2D0, pt2D0,
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PD0, pD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(P2D0, p2D0,
                           [](float px, float py, float pz) -> float { return RecoDecay::p2(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorD0, pVectorD0,
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(EtaD0, etaD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiD0, phiD0,
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(YD0, yD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::y(std::array{px, py, pz}, constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(ED0, eD0,
                           [](float px, float py, float pz) -> float { return RecoDecay::e(px, py, pz, constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(E2D0, e2D0,
                           [](float px, float py, float pz) -> float { return RecoDecay::e2(px, py, pz, constants::physics::MassD0); });
// secondary vertex
DECLARE_SOA_COLUMN(Chi2PCAD0, chi2PCAD0, float);
DECLARE_SOA_COLUMN(XSecondaryVertexD0, xSecondaryVertexD0, float);
DECLARE_SOA_COLUMN(YSecondaryVertexD0, ySecondaryVertexD0, float);
DECLARE_SOA_COLUMN(ZSecondaryVertexD0, zSecondaryVertexD0, float);
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertexD0, rSecondaryVertexD0,
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthD0, decayLengthD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYD0, decayLengthXYD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalisedD0, decayLengthNormalisedD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalisedD0, decayLengthXYNormalisedD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLengthD0, errorDecayLengthD0, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYD0, errorDecayLengthXYD0, float);
DECLARE_SOA_DYNAMIC_COLUMN(CPAD0, cpaD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CPAXYD0, cpaXYD0,
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(CtD0, ctD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), constants::physics::MassD0); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXYD0, impactParameterXYD0,
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });

// Columns only for D* properties
DECLARE_SOA_INDEX_COLUMN_FULL(ProngPi, prongPi, int, Tracks, ""); //! soft-pion index

// soft pion prong
DECLARE_SOA_COLUMN(ImpParamSoftPi, impParamSoftPi, float);
DECLARE_SOA_COLUMN(ImpParamZSoftPi, impParamZSoftPi, float);
DECLARE_SOA_COLUMN(ErrorImpParamSoftPi, errorImpParamSoftPi, float);
DECLARE_SOA_COLUMN(ErrorImpParamZSoftPi, errorImpParamZSoftPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(NormalisedImpParamSoftPi, normalisedImpParamSoftPi,
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_DYNAMIC_COLUMN(NormalisedImpParamZSoftPi, normalisedImpParamZSoftPi,
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxSoftPi, pxSoftPi, float);
DECLARE_SOA_COLUMN(PySoftPi, pySoftPi, float);
DECLARE_SOA_COLUMN(PzSoftPi, pzSoftPi, float);
DECLARE_SOA_COLUMN(SignSoftPi, signSoftPi, int8_t);
// Dstar momenta
DECLARE_SOA_EXPRESSION_COLUMN(PxDstar, pxDstar, float, 1.f * aod::hf_cand::pxProng0 + 1.f * aod::hf_cand::pxProng1 + 1.f * aod::hf_cand_dstar::pxSoftPi);
DECLARE_SOA_EXPRESSION_COLUMN(PyDstar, pyDstar, float, 1.f * aod::hf_cand::pyProng0 + 1.f * aod::hf_cand::pyProng1 + 1.f * aod::hf_cand_dstar::pySoftPi);
DECLARE_SOA_EXPRESSION_COLUMN(PzDstar, pzDstar, float, 1.f * aod::hf_cand::pzProng0 + 1.f * aod::hf_cand::pzProng1 + 1.f * aod::hf_cand_dstar::pzSoftPi);
// Inv Mass (accept mass array of size 3 {π , π, k})
DECLARE_SOA_DYNAMIC_COLUMN(InvMassDstar, invMassDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassPiPlus, constants::physics::MassKPlus}); });

DECLARE_SOA_DYNAMIC_COLUMN(InvMassAntiDstar, invMassAntiDstar,
                           [](float pxSoftPi, float pySoftPi, float pzSoftPi, float pxProng0, float pyProng0, float pzProng0, float pxProng1, float pyProng1, float pzProng1)
                             -> float { return RecoDecay::m(std::array{std::array{pxSoftPi, pySoftPi, pzSoftPi}, std::array{pxProng0, pyProng0, pzProng0}, std::array{pxProng1, pyProng1, pzProng1}}, std::array{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassPiPlus}); });

DECLARE_SOA_DYNAMIC_COLUMN(PtSoftPi, ptSoftPi, [](float pxSoftPi, float pySoftPi) -> float { return RecoDecay::pt(pxSoftPi, pySoftPi); });
DECLARE_SOA_DYNAMIC_COLUMN(PVecSoftPi, pVecSoftPi, [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });

// MC matching result:
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); //! reconstruction level
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t); //! generator level

DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); //! particle origin, generator level

enum DecayType {
  DstarToD0Pi = 0,
  D0ToPiK,
  NDstarDecayType
};

} // namespace hf_cand_dstar

/// D0 (table) from DStar
DECLARE_SOA_TABLE(HfD0FromDstarBase, "AOD", "HFD0FRMDSTR",
                  o2::soa::Index<>,
                  // gener columns
                  hf_cand::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0,
                  hf_cand_dstar::ErrorDecayLengthD0, hf_cand_dstar::ErrorDecayLengthXYD0,
                  hf_cand_dstar::Chi2PCAD0,
                  /* dynamic columns */ hf_cand_dstar::RSecondaryVertexD0<hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthXYD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0>,
                  hf_cand_dstar::DecayLengthNormalisedD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthD0>,
                  hf_cand_dstar::DecayLengthXYNormalisedD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthXYD0>,
                  /* prong 0 */ hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  /* prong 1 */ hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  /* prong 0 */ hf_cand::ImpactParameterZNormalised0<hf_cand::ImpactParameterZ0, hf_cand::ErrorImpactParameterZ0>,
                  /* prong 1 */ hf_cand::ImpactParameterZNormalised1<hf_cand::ImpactParameterZ1, hf_cand::ErrorImpactParameterZ1>,

                  // HFCAND_COLUMNS,
                  // 2-prong specific columns
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_dstar::ImpactParameterProductD0<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_dstar::ImpactParameterProngSqSumD0<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand_dstar::PtD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::Pt2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::PD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::P2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::PVectorD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::CPAD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::CPAXYD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::CtD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::ImpactParameterXYD0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ZSecondaryVertexD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::DeltaIPNormalisedMaxD0<collision::PosX, collision::PosY, hf_cand_dstar::XSecondaryVertexD0, hf_cand_dstar::YSecondaryVertexD0, hf_cand_dstar::ErrorDecayLengthXYD0, hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0, hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand_dstar::EtaD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::PhiD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0>,
                  hf_cand_dstar::YD0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::ED0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>,
                  hf_cand_dstar::E2D0<hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfD0FromDstarExt, HfD0FromDstarBase, "HFD0FRMDSTREXT",
                                hf_cand_dstar::PxD0, hf_cand_dstar::PyD0, hf_cand_dstar::PzD0);

using HfD0FromDstar = HfD0FromDstarExt;

DECLARE_SOA_TABLE(HfCandDstarBase, "AOD", "HFCANDDSTRBASE",
                  o2::soa::Index<>,
                  hf_cand::CollisionId,
                  // Primary vertex
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_dstar::ProngPiId,  // Index column to softPi track table
                  hf_track_index::ProngD0Id, // Index column points to Hf2Prongs table filled by indexSkimcreator
                  // Softpi
                  hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi,
                  hf_cand_dstar::SignSoftPi,
                  hf_cand_dstar::ImpParamSoftPi, hf_cand_dstar::ImpParamZSoftPi,
                  hf_cand_dstar::ErrorImpParamSoftPi, hf_cand_dstar::ErrorImpParamZSoftPi,
                  // Two pronges of D0
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  // Dynamic
                  hf_cand_dstar::PtSoftPi<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi>,
                  hf_cand_dstar::PVecSoftPi<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi>,
                  hf_cand_dstar::NormalisedImpParamSoftPi<hf_cand_dstar::ImpParamSoftPi, hf_cand_dstar::ErrorImpParamSoftPi>,
                  hf_cand_dstar::NormalisedImpParamZSoftPi<hf_cand_dstar::ImpParamZSoftPi, hf_cand_dstar::ErrorImpParamZSoftPi>,
                  hf_cand::Pt<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar>,
                  hf_cand::P<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::PVector<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::Eta<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::Phi<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar>,
                  hf_cand::Y<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::E<hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar>,
                  hf_cand::PtProng0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::Pt2Prong0<hf_cand::PxProng0, hf_cand::PyProng0>,
                  hf_cand::PVectorProng0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0>,
                  hf_cand::PtProng1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::Pt2Prong1<hf_cand::PxProng1, hf_cand::PyProng1>,
                  hf_cand::PVectorProng1<hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassDstar<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassAntiDstar<hf_cand_dstar::PxSoftPi, hf_cand_dstar::PySoftPi, hf_cand_dstar::PzSoftPi, hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMassD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMass2D0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::InvMass2D0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::CosThetaStarD0<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand_dstar::CosThetaStarD0Bar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandDstarExt, HfCandDstarBase, "HFCANDDSTREXT",
                                hf_cand_dstar::PxDstar, hf_cand_dstar::PyDstar, hf_cand_dstar::PzDstar);

using HfCandDstars = HfCandDstarExt;
using HfCandDstar = HfCandDstars::iterator;

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandDstarMcRec, "AOD", "HFCANDDSTRMCREC",
                  hf_cand_dstar::FlagMcMatchRec,
                  hf_cand_dstar::OriginMcRec,
                  hf_cand::PtBhadMotherPart,
                  hf_cand::PdgBhadMotherPart,
                  hf_cand::NTracksDecayed,
                  hf_cand::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandDstarMcGen, "AOD", "HFCANDDSTRMCGEN",
                  hf_cand_dstar::FlagMcMatchGen,
                  hf_cand_dstar::OriginMcGen,
                  hf_cand::IdxBhadMotherPart);

#undef HFCAND_COLUMNS

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLES_H_
