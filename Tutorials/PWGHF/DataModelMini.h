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

/// \file DataModelMini.h
/// \brief Mini version of the HF analysis chain
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef TUTORIALS_PWGHF_DATAMODELMINI_H_
#define TUTORIALS_PWGHF_DATAMODELMINI_H_

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <array>

namespace o2::aod
{
namespace hf_seltrack
{
// Track selection columns
DECLARE_SOA_COLUMN(IsSelProng, isSelProng, bool); //! prong selection flag
} // namespace hf_seltrack

// Track selection table
DECLARE_SOA_TABLE(HfTSelTrack, "AOD", "HFTSELTRACK", //! track selection table
                  hf_seltrack::IsSelProng);

namespace hf_track_index
{
// Track index skim columns
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, Tracks, "_0"); //! prong 0
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, Tracks, "_1"); //! prong 1
} // namespace hf_track_index

// Track index skim table
DECLARE_SOA_TABLE(HfT2Prongs, "AOD", "HFT2PRONG", //! table with prong indices
                  o2::soa::Index<>,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);

// 2-prong decay properties
namespace hf_cand_prong2
{
// Candidate columns
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! collision
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float); //! x coordinate of the secondary vertex [cm]
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float); //! y coordinate of the secondary vertex [cm]
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float); //! z coordinate of the secondary vertex [cm]
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, //! radius of the secondary vertex [cm]
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //! px of prong 0 [GeV/c]
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //! py of prong 0 [GeV/c]
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //! pz of prong 0 [GeV/c]
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //! pt of prong 0 [GeV/c]
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //! px of prong 1 [GeV/c]
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //! py of prong 1 [GeV/c]
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //! pz of prong 1 [GeV/c]
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //! pt of prong 1 [GeV/c]
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //! decay length of candidate [cm]
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! pt of candidate [GeV/c]
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! px of candidate [GeV/c]
                              float, 1.f * pxProng0 + 1.f * pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! py of candidate [GeV/c]
                              float, 1.f * pyProng0 + 1.f * pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! pz of candidate [GeV/c]
                              float, 1.f * pzProng0 + 1.f * pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! invariant mass of candidate [GeV/c^2]
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(Cpa, cpa, //! cosine of pointing angle of candidate
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
} // namespace hf_cand_prong2

// Candidate table
DECLARE_SOA_TABLE(HfTCand2ProngBase, "AOD", "HFTCAND2PBASE", //! 2-prong candidate table
                  o2::soa::Index<>,
                  hf_cand_prong2::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex,
                  /* dynamic columns */ hf_cand_prong2::RSecondaryVertex<hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex>,
                  hf_cand_prong2::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex>,
                  /* prong 0 */ hf_cand_prong2::PtProng0<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0>,
                  hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0,
                  /* prong 1 */ hf_cand_prong2::PtProng1<hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1>,
                  hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  /* dynamic columns */
                  hf_cand_prong2::M<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0, hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand_prong2::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex, hf_cand_prong2::Px, hf_cand_prong2::Py, hf_cand_prong2::Pz>,
                  hf_cand_prong2::Pt<hf_cand_prong2::Px, hf_cand_prong2::Py>);

// Extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfTCand2ProngExt, HfTCand2ProngBase, "HFTCAND2PEXT", //! extension table for the 2-prong candidate table
                                hf_cand_prong2::Px, hf_cand_prong2::Py, hf_cand_prong2::Pz);

using HfTCand2Prong = HfTCand2ProngExt;

namespace hf_selcandidate_d0
{
// Candidate selection columns
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);       //! selection flag for D0
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int); //! selection flag for D0bar
} // namespace hf_selcandidate_d0

// Candidate selection table
DECLARE_SOA_TABLE(HfTSelD0, "AOD", "HFTSELD0", //! table with D0 selection flags
                  hf_selcandidate_d0::IsSelD0,
                  hf_selcandidate_d0::IsSelD0bar);
} // namespace o2::aod

#endif // TUTORIALS_PWGHF_DATAMODELMINI_H_
