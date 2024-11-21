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

/// \file DerivedTables.h
/// \brief Definitions of derived tables produced by derived-data creators
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_DATAMODEL_DERIVEDTABLES_H_
#define PWGHF_DATAMODEL_DERIVEDTABLES_H_

#include <vector>

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

namespace o2::aod
{
// basic species:
// D0 -> K- + pi+ (done)
// Lc -> pi+ K- p (existing 3P table to be renamed Lc)
// D+ -> K- + pi+ + pi+ (3P table with adapted PID columns)
// Ds+ -> K- + K+ + pi+ (3P table with adapted PID columns)
// composite species
// B0 -> D- + pi+
// B+ -> D0 + pi+ (in progress)
// D*+ -> D0 + pi+
constexpr uint MarkerBase = 2;
constexpr uint MarkerD0 = 3;
constexpr uint Marker3P = 4;
constexpr uint MarkerBplus = 5;
constexpr uint MarkerB0 = 6;

// ================
// Collision tables
// ================

// Basic collision properties
namespace hf_coll_base
{
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int8_t);      //! collision rejection flag
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);                 //! FT0M multiplicity
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);                 //! FT0A centrality percentile
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);                 //! FT0C centrality percentile
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);                 //! FT0M centrality percentile
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);                 //! FV0A centrality percentile
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);                 //! FDDM centrality percentile
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float); //! z-equalised barrel multiplicity
} // namespace hf_coll_base

// base

DECLARE_SOA_TABLE_STAGED(HfCollBases, "HFCOLLBASE", //! Table with basic collision info
                         o2::soa::Index<>,
                         collision::PosX,
                         collision::PosY,
                         collision::PosZ,
                         collision::NumContrib,
                         hf_coll_base::CentFT0A,
                         hf_coll_base::CentFT0C,
                         hf_coll_base::CentFT0M,
                         hf_coll_base::CentFV0A,
                         hf_coll_base::MultZeqNTracksPV,
                         // hf_coll_base::IsEventReject,
                         // bc::RunNumber,
                         o2::soa::Marker<MarkerBase>);

using HfCollBase = HfCollBases::iterator;

DECLARE_SOA_TABLE_STAGED(HfCollIds, "HFCOLLID", //! Table with original global indices of collisions
                         hf_cand::CollisionId,
                         o2::soa::Marker<MarkerBase>);

// D0 (to be replaced by base version)

DECLARE_SOA_TABLE_STAGED(HfD0CollBases, "HFD0COLLBASE", //! Table with basic collision info
                         o2::soa::Index<>,
                         collision::PosX,
                         collision::PosY,
                         collision::PosZ,
                         collision::NumContrib,
                         hf_coll_base::CentFT0A,
                         hf_coll_base::CentFT0C,
                         hf_coll_base::CentFT0M,
                         hf_coll_base::CentFV0A,
                         hf_coll_base::MultZeqNTracksPV,
                         // hf_coll_base::IsEventReject,
                         // bc::RunNumber,
                         o2::soa::Marker<MarkerD0>);

using HfD0CollBase = HfD0CollBases::iterator;
using StoredHfD0CollBase = StoredHfD0CollBases::iterator;

DECLARE_SOA_TABLE_STAGED(HfD0CollIds, "HFD0COLLID", //! Table with original global indices of collisions
                         hf_cand::CollisionId,
                         o2::soa::Marker<MarkerD0>);

// 3-prong decays (to be replaced by base version)

DECLARE_SOA_TABLE_STAGED(Hf3PCollBases, "HF3PCOLLBASE", //! Table with basic collision info
                         o2::soa::Index<>,
                         collision::PosX,
                         collision::PosY,
                         collision::PosZ,
                         collision::NumContrib,
                         hf_coll_base::CentFT0A,
                         hf_coll_base::CentFT0C,
                         hf_coll_base::CentFT0M,
                         hf_coll_base::CentFV0A,
                         hf_coll_base::MultZeqNTracksPV,
                         // hf_coll_base::IsEventReject,
                         // bc::RunNumber,
                         o2::soa::Marker<Marker3P>);

using Hf3PCollBase = Hf3PCollBases::iterator;
using StoredHf3PCollBase = StoredHf3PCollBases::iterator;

DECLARE_SOA_TABLE_STAGED(Hf3PCollIds, "HF3PCOLLID", //! Table with original global indices of collisions
                         hf_cand::CollisionId,
                         o2::soa::Marker<Marker3P>);

// ===================
// MC collision tables
// ===================

// MC collision columns
namespace hf_mc_coll
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);      //! original global index of the MC collision
DECLARE_SOA_ARRAY_INDEX_COLUMN(HfCollBase, hfCollBases); //! collision index array pointing to the derived reconstructed collisions for D0 candidates
namespace der_d0
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(HfD0CollBase, hfCollBases); //! collision index array pointing to the derived reconstructed collisions for D0 candidates
}
namespace der_3p
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(Hf3PCollBase, hfCollBases); //! collision index array pointing to the derived reconstructed collisions for 3-prong candidates
}
} // namespace hf_mc_coll

// base

DECLARE_SOA_TABLE_STAGED(HfMcCollBases, "HFMCCOLLBASE", //! Table with basic MC collision info
                         o2::soa::Index<>,
                         mccollision::PosX,
                         mccollision::PosY,
                         mccollision::PosZ,
                         o2::soa::Marker<MarkerBase>);

using HfMcCollBase = HfMcCollBases::iterator;

DECLARE_SOA_TABLE_STAGED(HfMcCollIds, "HFMCCOLLID", //! Table with original global indices of MC collisions
                         hf_mc_coll::McCollisionId,
                         o2::soa::Marker<MarkerBase>);

DECLARE_SOA_TABLE_STAGED(HfMcRCollIds, "HFMCRCOLLID", //! Table with indices pointing to the derived reconstructed-collision table
                         hf_mc_coll::HfCollBaseIds);

// D0

DECLARE_SOA_TABLE_STAGED(HfD0McCollBases, "HFD0MCCOLLBASE", //! Table with basic MC collision info
                         o2::soa::Index<>,
                         mccollision::PosX,
                         mccollision::PosY,
                         mccollision::PosZ,
                         o2::soa::Marker<MarkerD0>);

using HfD0McCollBase = HfD0McCollBases::iterator;
using StoredHfD0McCollBase = StoredHfD0McCollBases::iterator;

DECLARE_SOA_TABLE_STAGED(HfD0McCollIds, "HFD0MCCOLLID", //! Table with original global indices of MC collisions
                         hf_mc_coll::McCollisionId,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0McRCollIds, "HFD0MCRCOLLID", //! Table with indices pointing to the derived reconstructed-collision table
                         hf_mc_coll::der_d0::HfD0CollBaseIds);

// 3-prong decays

DECLARE_SOA_TABLE_STAGED(Hf3PMcCollBases, "HF3PMCCOLLBASE", //! Table with basic MC collision info
                         o2::soa::Index<>,
                         mccollision::PosX,
                         mccollision::PosY,
                         mccollision::PosZ,
                         o2::soa::Marker<Marker3P>);

using Hf3PMcCollBase = Hf3PMcCollBases::iterator;
using StoredHf3PMcCollBase = StoredHf3PMcCollBases::iterator;

DECLARE_SOA_TABLE_STAGED(Hf3PMcCollIds, "HF3PMCCOLLID", //! Table with original global indices of MC collisions
                         hf_mc_coll::McCollisionId,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PMcRCollIds, "HF3PMCRCOLLID", //! Table with indices pointing to the derived reconstructed-collision table
                         hf_mc_coll::der_3p::Hf3PCollBaseIds);

// ================
// Candidate tables
// ================

// Basic candidate properties
namespace hf_cand_base
{
namespace der_d0
{
DECLARE_SOA_INDEX_COLUMN(HfD0CollBase, hfCollBase); //! collision index pointing to the derived collision table for D0 candidates
}
namespace der_bplus
{
DECLARE_SOA_INDEX_COLUMN(HfCollBase, hfCollBase); //! collision index pointing to the derived collision table for B+ candidates
}
namespace der_3p
{
DECLARE_SOA_INDEX_COLUMN(Hf3PCollBase, hfCollBase); //! collision index pointing to the derived collision table for 3-prong candidates
}
DECLARE_SOA_COLUMN(Eta, eta, float); //! pseudorapidity
DECLARE_SOA_COLUMN(M, m, float);     //! invariant mass
DECLARE_SOA_COLUMN(Phi, phi, float); //! azimuth
DECLARE_SOA_COLUMN(Pt, pt, float);   //! transverse momentum
DECLARE_SOA_COLUMN(Y, y, float);     //! rapidity
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,   //! px
                           [](float pt, float phi) -> float { return RecoDecayPtEtaPhi::px(pt, phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! py
                           [](float pt, float phi) -> float { return RecoDecayPtEtaPhi::py(pt, phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! px
                           [](float pt, float eta) -> float { return RecoDecayPtEtaPhi::pz(pt, eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! momentum
                           [](float pt, float eta) -> float { return RecoDecayPtEtaPhi::p(pt, eta); });
} // namespace hf_cand_base

// Candidate properties used for selection
namespace hf_cand_par
{
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);                             //! cosine of theta star
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! cosine of pointing angle
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! cosine of pointing angle in the transverse plane
DECLARE_SOA_COLUMN(Ct, ct, float);                                                 //! proper lifetime times c
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! decay length
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! decay length divided by its uncertainty
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! decay length in the transverse plane
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! decay length in the transverse plane divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float); //! impact parameter of prong 0 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float); //! impact parameter of prong 1 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float); //! impact parameter of prong 2 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);         //! product of impact parameters of prong 0 and prong 1
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);             //! see RecoDecay::maxNormalisedDeltaIP
DECLARE_SOA_COLUMN(PProng0, pProng0, float);                                       //! momentum magnitude of prong 0
DECLARE_SOA_COLUMN(PProng1, pProng1, float);                                       //! momentum magnitude of prong 1
DECLARE_SOA_COLUMN(PProng2, pProng2, float);                                       //! momentum magnitude of prong 2
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                                     //! transverse momentum of prong 0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                                     //! transverse momentum of prong 1
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);                                     //! transverse momentum of prong 2
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);                     //! distance of the secondary vertex from the z axis
// TOF
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);
DECLARE_SOA_COLUMN(NSigTofKaExpPi, nSigTofKaExpPi, float);
DECLARE_SOA_COLUMN(NSigTofKaExpKa, nSigTofKaExpKa, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);
DECLARE_SOA_COLUMN(NSigTofPiExpPi, nSigTofPiExpPi, float);
DECLARE_SOA_COLUMN(NSigTofPiExpKa, nSigTofPiExpKa, float);
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);
DECLARE_SOA_COLUMN(NSigTofPr1, nSigTofPr1, float);
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);
// TPC
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);
DECLARE_SOA_COLUMN(NSigTpcKaExpPi, nSigTpcKaExpPi, float);
DECLARE_SOA_COLUMN(NSigTpcKaExpKa, nSigTpcKaExpKa, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);
DECLARE_SOA_COLUMN(NSigTpcPiExpPi, nSigTpcPiExpPi, float);
DECLARE_SOA_COLUMN(NSigTpcPiExpKa, nSigTpcPiExpKa, float);
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);
DECLARE_SOA_COLUMN(NSigTpcPr1, nSigTpcPr1, float);
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);
// TPC+TOF
DECLARE_SOA_COLUMN(NSigTpcTofKa0, nSigTpcTofKa0, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa2, nSigTpcTofKa2, float);
DECLARE_SOA_COLUMN(NSigTpcTofKaExpPi, nSigTpcTofKaExpPi, float);
DECLARE_SOA_COLUMN(NSigTpcTofKaExpKa, nSigTpcTofKaExpKa, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPiExpPi, nSigTpcTofPiExpPi, float);
DECLARE_SOA_COLUMN(NSigTpcTofPiExpKa, nSigTpcTofPiExpKa, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr0, nSigTpcTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr1, nSigTpcTofPr1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr2, nSigTpcTofPr2, float);
} // namespace hf_cand_par

// Candidate properties of the charm daughter candidate used for selection of the beauty candidate
// Copy of hf_cand_par with "Charm" suffix to make it joinable with the beauty candidate table.
// We don't want to link the charm candidate table because we want to avoid producing it.
namespace hf_cand_par_charm
{
DECLARE_SOA_COLUMN(CosThetaStarCharm, cosThetaStarCharm, float);                             //! cosine of theta star
DECLARE_SOA_COLUMN(CpaCharm, cpaCharm, float);                                               //! cosine of pointing angle
DECLARE_SOA_COLUMN(CpaXYCharm, cpaXYCharm, float);                                           //! cosine of pointing angle in the transverse plane
DECLARE_SOA_COLUMN(CtCharm, ctCharm, float);                                                 //! proper lifetime times c
DECLARE_SOA_COLUMN(DecayLengthCharm, decayLengthCharm, float);                               //! decay length
DECLARE_SOA_COLUMN(DecayLengthNormalisedCharm, decayLengthNormalisedCharm, float);           //! decay length divided by its uncertainty
DECLARE_SOA_COLUMN(DecayLengthXYCharm, decayLengthXYCharm, float);                           //! decay length in the transverse plane
DECLARE_SOA_COLUMN(DecayLengthXYNormalisedCharm, decayLengthXYNormalisedCharm, float);       //! decay length in the transverse plane divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameter0Charm, impactParameter0Charm, float);                     //! impact parameter of prong 0
DECLARE_SOA_COLUMN(ImpactParameter1Charm, impactParameter1Charm, float);                     //! impact parameter of prong 1
DECLARE_SOA_COLUMN(ImpactParameterNormalised0Charm, impactParameterNormalised0Charm, float); //! impact parameter of prong 0 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterNormalised1Charm, impactParameterNormalised1Charm, float); //! impact parameter of prong 1 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterNormalised2Charm, impactParameterNormalised2Charm, float); //! impact parameter of prong 2 divided by its uncertainty
DECLARE_SOA_COLUMN(ImpactParameterProductCharm, impactParameterProductCharm, float);         //! product of impact parameters of prong 0 and prong 1
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIPCharm, maxNormalisedDeltaIPCharm, float);             //! see RecoDecay::maxNormalisedDeltaIP
DECLARE_SOA_COLUMN(PProng0Charm, pProng0Charm, float);                                       //! momentum magnitude of prong 0
DECLARE_SOA_COLUMN(PProng1Charm, pProng1Charm, float);                                       //! momentum magnitude of prong 1
DECLARE_SOA_COLUMN(PProng2Charm, pProng2Charm, float);                                       //! momentum magnitude of prong 2
DECLARE_SOA_COLUMN(PtProng0Charm, ptProng0Charm, float);                                     //! transverse momentum of prong 0
DECLARE_SOA_COLUMN(PtProng1Charm, ptProng1Charm, float);                                     //! transverse momentum of prong 1
DECLARE_SOA_COLUMN(PtProng2Charm, ptProng2Charm, float);                                     //! transverse momentum of prong 2
DECLARE_SOA_COLUMN(RSecondaryVertexCharm, rSecondaryVertexCharm, float);                     //! distance of the secondary vertex from the z axis
// TOF
DECLARE_SOA_COLUMN(NSigTofKa0Charm, nSigTofKa0Charm, float);
DECLARE_SOA_COLUMN(NSigTofKa1Charm, nSigTofKa1Charm, float);
DECLARE_SOA_COLUMN(NSigTofKa2Charm, nSigTofKa2Charm, float);
DECLARE_SOA_COLUMN(NSigTofKaExpPiCharm, nSigTofKaExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTofKaExpKaCharm, nSigTofKaExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTofPi0Charm, nSigTofPi0Charm, float);
DECLARE_SOA_COLUMN(NSigTofPi1Charm, nSigTofPi1Charm, float);
DECLARE_SOA_COLUMN(NSigTofPi2Charm, nSigTofPi2Charm, float);
DECLARE_SOA_COLUMN(NSigTofPiExpPiCharm, nSigTofPiExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTofPiExpKaCharm, nSigTofPiExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTofPr0Charm, nSigTofPr0Charm, float);
DECLARE_SOA_COLUMN(NSigTofPr1Charm, nSigTofPr1Charm, float);
DECLARE_SOA_COLUMN(NSigTofPr2Charm, nSigTofPr2Charm, float);
// TPC
DECLARE_SOA_COLUMN(NSigTpcKa0Charm, nSigTpcKa0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcKa1Charm, nSigTpcKa1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcKa2Charm, nSigTpcKa2Charm, float);
DECLARE_SOA_COLUMN(NSigTpcKaExpPiCharm, nSigTpcKaExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTpcKaExpKaCharm, nSigTpcKaExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTpcPi0Charm, nSigTpcPi0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcPi1Charm, nSigTpcPi1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcPi2Charm, nSigTpcPi2Charm, float);
DECLARE_SOA_COLUMN(NSigTpcPiExpPiCharm, nSigTpcPiExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTpcPiExpKaCharm, nSigTpcPiExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTpcPr0Charm, nSigTpcPr0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcPr1Charm, nSigTpcPr1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcPr2Charm, nSigTpcPr2Charm, float);
// TPC+TOF
DECLARE_SOA_COLUMN(NSigTpcTofKa0Charm, nSigTpcTofKa0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1Charm, nSigTpcTofKa1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa2Charm, nSigTpcTofKa2Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofKaExpPiCharm, nSigTpcTofKaExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTpcTofKaExpKaCharm, nSigTpcTofKaExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0Charm, nSigTpcTofPi0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi1Charm, nSigTpcTofPi1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi2Charm, nSigTpcTofPi2Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPiExpPiCharm, nSigTpcTofPiExpPiCharm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPiExpKaCharm, nSigTpcTofPiExpKaCharm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr0Charm, nSigTpcTofPr0Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr1Charm, nSigTpcTofPr1Charm, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr2Charm, nSigTpcTofPr2Charm, float);
} // namespace hf_cand_par_charm

// Candidate selection flags
namespace hf_cand_sel
{
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t); //! bitmap of the selected candidate type
}

// Candidate MC columns
namespace hf_cand_mc
{
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         //! flag for reconstruction level matching
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t); //! swapping of the prongs order
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                  //! ML score for signal class
DECLARE_SOA_COLUMN(MlScoreBkg, mlScoreBkg, float);                  //! ML score for background class
DECLARE_SOA_COLUMN(MlScorePrompt, mlScorePrompt, float);            //! ML score for prompt class
DECLARE_SOA_COLUMN(MlScoreNonPrompt, mlScoreNonPrompt, float);      //! ML score for non-prompt class
DECLARE_SOA_COLUMN(MlScores, mlScores, std::vector<float>);         //! vector of ML scores
} // namespace hf_cand_mc

// Candidate MC columns of the charm daughter
namespace hf_cand_mc_charm
{
DECLARE_SOA_COLUMN(FlagMcMatchRecCharm, flagMcMatchRecCharm, int8_t);         //! flag for reconstruction level matching
DECLARE_SOA_COLUMN(OriginMcRecCharm, originMcRecCharm, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(IsCandidateSwappedCharm, isCandidateSwappedCharm, int8_t); //! swapping of the prongs order
DECLARE_SOA_COLUMN(FlagMcDecayChanRecCharm, flagMcDecayChanRecCharm, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(MlScoreSigCharm, mlScoreSigCharm, float);                  //! ML score for signal class
DECLARE_SOA_COLUMN(MlScoreBkgCharm, mlScoreBkgCharm, float);                  //! ML score for background class
DECLARE_SOA_COLUMN(MlScorePromptCharm, mlScorePromptCharm, float);            //! ML score for prompt class
DECLARE_SOA_COLUMN(MlScoreNonPromptCharm, mlScoreNonPromptCharm, float);      //! ML score for non-prompt class
DECLARE_SOA_COLUMN(MlScoresCharm, mlScoresCharm, std::vector<float>);         //! vector of ML scores
} // namespace hf_cand_mc_charm

// D0

DECLARE_SOA_TABLE_STAGED(HfD0Bases, "HFD0BASE", //! Table with basic candidate properties used in the analyses
                         o2::soa::Index<>,
                         hf_cand_base::der_d0::HfD0CollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::M,
                         hf_cand_base::Y,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<MarkerD0>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE_STAGED(HfD0Pars, "HFD0PAR", //! Table with candidate properties used for selection
                         hf_cand::Chi2PCA,
                         hf_cand_par::Cpa,
                         hf_cand_par::CpaXY,
                         hf_cand_par::DecayLength,
                         hf_cand_par::DecayLengthXY,
                         hf_cand_par::DecayLengthNormalised,
                         hf_cand_par::DecayLengthXYNormalised,
                         hf_cand_par::PtProng0,
                         hf_cand_par::PtProng1,
                         hf_cand::ImpactParameter0,
                         hf_cand::ImpactParameter1,
                         hf_cand_par::ImpactParameterNormalised0,
                         hf_cand_par::ImpactParameterNormalised1,
                         hf_cand_par::NSigTpcPiExpPi,
                         hf_cand_par::NSigTofPiExpPi,
                         hf_cand_par::NSigTpcTofPiExpPi,
                         hf_cand_par::NSigTpcKaExpPi,
                         hf_cand_par::NSigTofKaExpPi,
                         hf_cand_par::NSigTpcTofKaExpPi,
                         hf_cand_par::NSigTpcPiExpKa,
                         hf_cand_par::NSigTofPiExpKa,
                         hf_cand_par::NSigTpcTofPiExpKa,
                         hf_cand_par::NSigTpcKaExpKa,
                         hf_cand_par::NSigTofKaExpKa,
                         hf_cand_par::NSigTpcTofKaExpKa,
                         hf_cand_par::MaxNormalisedDeltaIP,
                         hf_cand_par::ImpactParameterProduct,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0ParEs, "HFD0PARE", //! Table with additional candidate properties used for selection
                         hf_cand::XSecondaryVertex,
                         hf_cand::YSecondaryVertex,
                         hf_cand::ZSecondaryVertex,
                         hf_cand::ErrorDecayLength,
                         hf_cand::ErrorDecayLengthXY,
                         hf_cand::KfTopolChi2OverNdf,
                         hf_cand_par::RSecondaryVertex,
                         hf_cand_par::PProng0,
                         hf_cand_par::PProng1,
                         hf_cand::PxProng0,
                         hf_cand::PyProng0,
                         hf_cand::PzProng0,
                         hf_cand::PxProng1,
                         hf_cand::PyProng1,
                         hf_cand::PzProng1,
                         hf_cand::ErrorImpactParameter0,
                         hf_cand::ErrorImpactParameter1,
                         hf_cand_par::CosThetaStar,
                         hf_cand_par::Ct,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0Sels, "HFD0SEL", //! Table with candidate selection flags
                         hf_cand_sel::CandidateSelFlag,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0Mls, "HFD0ML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScores,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0Ids, "HFD0ID", //! Table with original global indices for candidates
                         hf_cand::CollisionId,
                         hf_track_index::Prong0Id,
                         hf_track_index::Prong1Id,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0Mcs, "HFD0MC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         o2::soa::Marker<MarkerD0>);

// B+

DECLARE_SOA_TABLE_STAGED(HfBplusBases, "HFBPBASE", //! Table with basic candidate properties used in the analyses
                         o2::soa::Index<>,
                         hf_cand_base::der_bplus::HfCollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::M,
                         hf_cand_base::Y,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<MarkerBplus>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE_STAGED(HfBplusPars, "HFBPPAR", //! Table with candidate properties used for selection
                         hf_cand::Chi2PCA,
                         hf_cand_par::Cpa,
                         hf_cand_par::CpaXY,
                         hf_cand_par::DecayLength,
                         hf_cand_par::DecayLengthXY,
                         hf_cand_par::DecayLengthNormalised,
                         hf_cand_par::DecayLengthXYNormalised,
                         hf_cand_par::PtProng0,
                         hf_cand_par::PtProng1,
                         hf_cand::ImpactParameter0,
                         hf_cand::ImpactParameter1,
                         hf_cand_par::ImpactParameterNormalised0,
                         hf_cand_par::ImpactParameterNormalised1,
                         hf_cand_par::NSigTpcPiExpPi,
                         hf_cand_par::NSigTofPiExpPi,
                         hf_cand_par::NSigTpcTofPiExpPi,
                         hf_cand_par::NSigTpcKaExpPi,
                         hf_cand_par::NSigTofKaExpPi,
                         hf_cand_par::NSigTpcTofKaExpPi,
                         hf_cand_par::MaxNormalisedDeltaIP,
                         hf_cand_par::ImpactParameterProduct,
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusParD0s, "HFBPPARD0", //! Table with D0 candidate properties used for selection of B+
                         hf_cand_par_charm::CpaCharm,
                         hf_cand_par_charm::DecayLengthCharm,
                         hf_cand_par_charm::ImpactParameter0Charm,
                         hf_cand_par_charm::ImpactParameter1Charm,
                         hf_cand_par_charm::ImpactParameterProductCharm,
                         hf_cand_par_charm::NSigTpcPiExpPiCharm,
                         hf_cand_par_charm::NSigTofPiExpPiCharm,
                         hf_cand_par_charm::NSigTpcTofPiExpPiCharm,
                         hf_cand_par_charm::NSigTpcKaExpPiCharm,
                         hf_cand_par_charm::NSigTofKaExpPiCharm,
                         hf_cand_par_charm::NSigTpcTofKaExpPiCharm,
                         hf_cand_par_charm::NSigTpcPiExpKaCharm,
                         hf_cand_par_charm::NSigTofPiExpKaCharm,
                         hf_cand_par_charm::NSigTpcTofPiExpKaCharm,
                         hf_cand_par_charm::NSigTpcKaExpKaCharm,
                         hf_cand_par_charm::NSigTofKaExpKaCharm,
                         hf_cand_par_charm::NSigTpcTofKaExpKaCharm);

DECLARE_SOA_TABLE_STAGED(HfBplusParEs, "HFBPPARE", //! Table with additional candidate properties used for selection
                         hf_cand::XSecondaryVertex,
                         hf_cand::YSecondaryVertex,
                         hf_cand::ZSecondaryVertex,
                         hf_cand::ErrorDecayLength,
                         hf_cand::ErrorDecayLengthXY,
                         hf_cand_par::RSecondaryVertex,
                         hf_cand_par::PProng1,
                         hf_cand::PxProng1,
                         hf_cand::PyProng1,
                         hf_cand::PzProng1,
                         hf_cand::ErrorImpactParameter1,
                         hf_cand_par::CosThetaStar,
                         hf_cand_par::Ct,
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusMls, "HFBPML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScoreSig,
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusMlD0s, "HFBPMLD0", //! Table with D0 candidate selection ML scores
                         hf_cand_mc_charm::MlScoresCharm,
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusIds, "HFBPID", //! Table with original global indices for candidates
                         hf_cand::CollisionId,
                         hf_track_index::Prong0Id, // D0 prong 0
                         hf_track_index::Prong1Id, // D0 prong 1
                         hf_track_index::Prong2Id, // bachelor pion
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusMcs, "HFBPMC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         o2::soa::Marker<MarkerBplus>);

// 3-prong decays

DECLARE_SOA_TABLE_STAGED(Hf3PBases, "HF3PBASE", //! Table with basic candidate properties used in the analyses
                         o2::soa::Index<>,
                         hf_cand_base::der_3p::Hf3PCollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::M,
                         hf_cand_base::Y,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<Marker3P>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE_STAGED(Hf3PPars, "HF3PPAR", //! Table with candidate properties used for selection
                         hf_cand::Chi2PCA,
                         hf_cand::NProngsContributorsPV,
                         hf_cand_par::Cpa,
                         hf_cand_par::CpaXY,
                         hf_cand_par::DecayLength,
                         hf_cand_par::DecayLengthXY,
                         hf_cand_par::DecayLengthNormalised,
                         hf_cand_par::DecayLengthXYNormalised,
                         hf_cand_par::PtProng0,
                         hf_cand_par::PtProng1,
                         hf_cand_par::PtProng2,
                         hf_cand::ImpactParameter0,
                         hf_cand::ImpactParameter1,
                         hf_cand::ImpactParameter2,
                         hf_cand_par::ImpactParameterNormalised0,
                         hf_cand_par::ImpactParameterNormalised1,
                         hf_cand_par::ImpactParameterNormalised2,
                         hf_cand_par::NSigTpcPi0,
                         hf_cand_par::NSigTpcPr0,
                         hf_cand_par::NSigTofPi0,
                         hf_cand_par::NSigTofPr0,
                         hf_cand_par::NSigTpcTofPi0,
                         hf_cand_par::NSigTpcTofPr0,
                         hf_cand_par::NSigTpcKa1,
                         hf_cand_par::NSigTofKa1,
                         hf_cand_par::NSigTpcTofKa1,
                         hf_cand_par::NSigTpcPi2,
                         hf_cand_par::NSigTpcPr2,
                         hf_cand_par::NSigTofPi2,
                         hf_cand_par::NSigTofPr2,
                         hf_cand_par::NSigTpcTofPi2,
                         hf_cand_par::NSigTpcTofPr2,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PParEs, "HF3PPARE", //! Table with additional candidate properties used for selection
                         hf_cand::XSecondaryVertex,
                         hf_cand::YSecondaryVertex,
                         hf_cand::ZSecondaryVertex,
                         hf_cand::ErrorDecayLength,
                         hf_cand::ErrorDecayLengthXY,
                         hf_cand_par::RSecondaryVertex,
                         hf_cand_par::PProng0,
                         hf_cand_par::PProng1,
                         hf_cand_par::PProng2,
                         hf_cand::PxProng0,
                         hf_cand::PyProng0,
                         hf_cand::PzProng0,
                         hf_cand::PxProng1,
                         hf_cand::PyProng1,
                         hf_cand::PzProng1,
                         hf_cand::PxProng2,
                         hf_cand::PyProng2,
                         hf_cand::PzProng2,
                         hf_cand::ErrorImpactParameter0,
                         hf_cand::ErrorImpactParameter1,
                         hf_cand::ErrorImpactParameter2,
                         hf_cand_par::Ct,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PSels, "HF3PSEL", //! Table with candidate selection flags
                         hf_cand_sel::CandidateSelFlag,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PMls, "HF3PML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScores,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PIds, "HF3PID", //! Table with original global indices for candidates
                         hf_cand::CollisionId,
                         hf_track_index::Prong0Id,
                         hf_track_index::Prong1Id,
                         hf_track_index::Prong2Id,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PMcs, "HF3PMC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         hf_cand_mc::IsCandidateSwapped,
                         o2::soa::Marker<Marker3P>);

// ==================
// MC particle tables
// ==================

// MC particle columns
namespace hf_mc_particle
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision); //! MC collision of this particle
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);   //! MC particle
namespace der_d0
{
DECLARE_SOA_INDEX_COLUMN(HfD0McCollBase, hfMcCollBase); //! collision index pointing to the derived MC collision table for D0 candidates
}
namespace der_bplus
{
DECLARE_SOA_INDEX_COLUMN(HfMcCollBase, hfMcCollBase); //! collision index pointing to the derived MC collision table for B+ candidates
}
namespace der_3p
{
DECLARE_SOA_INDEX_COLUMN(Hf3PMcCollBase, hfMcCollBase); //! collision index pointing to the derived MC collision table for 3-prong candidates
}
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! flag for generator level matching
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, generator level
} // namespace hf_mc_particle

// D0

DECLARE_SOA_TABLE_STAGED(HfD0PBases, "HFD0PBASE", //! Table with MC particle info
                         o2::soa::Index<>,
                         hf_mc_particle::der_d0::HfD0McCollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::Y,
                         hf_mc_particle::FlagMcMatchGen,
                         hf_mc_particle::OriginMcGen,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0PIds, "HFD0PID", //! Table with original global indices for MC particles
                         hf_mc_particle::McCollisionId,
                         hf_mc_particle::McParticleId,
                         o2::soa::Marker<MarkerD0>);

// B+

DECLARE_SOA_TABLE_STAGED(HfBplusPBases, "HFBPPBASE", //! Table with MC particle info
                         o2::soa::Index<>,
                         hf_mc_particle::der_bplus::HfMcCollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::Y,
                         hf_mc_particle::FlagMcMatchGen,
                         hf_mc_particle::OriginMcGen,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<MarkerBplus>);

DECLARE_SOA_TABLE_STAGED(HfBplusPIds, "HFBPPID", //! Table with original global indices for MC particles
                         hf_mc_particle::McCollisionId,
                         hf_mc_particle::McParticleId,
                         o2::soa::Marker<MarkerBplus>);

// 3-prong decays

DECLARE_SOA_TABLE_STAGED(Hf3PPBases, "HF3PPBASE", //! Table with MC particle info
                         o2::soa::Index<>,
                         hf_mc_particle::der_3p::Hf3PMcCollBaseId,
                         hf_cand_base::Pt,
                         hf_cand_base::Eta,
                         hf_cand_base::Phi,
                         hf_cand_base::Y,
                         hf_mc_particle::FlagMcMatchGen,
                         hf_mc_particle::OriginMcGen,
                         hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                         hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                         hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                         o2::soa::Marker<Marker3P>);

DECLARE_SOA_TABLE_STAGED(Hf3PPIds, "HF3PPID", //! Table with original global indices for MC particles
                         hf_mc_particle::McCollisionId,
                         hf_mc_particle::McParticleId,
                         o2::soa::Marker<Marker3P>);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLES_H_
