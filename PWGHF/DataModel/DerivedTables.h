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

#include "PWGLF/DataModel/mcCentrality.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

namespace o2::aod
{
// basic species:
// D0 → K− π+
// Λc → p K− π+
// D+ → K− π+ π+
// Ds+ → K− K+ π+ (todo)

// composite species
// B0 → D− π+ (todo)
// B+ → D0 π+
// D*+ → D0 π+ (todo)

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

namespace hf_mc_coll
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision); //! original global index of the MC collision
} // namespace hf_mc_coll

// Declares the base table with reconstructed collisions (CollBases) and joinable tables (CollIds).
#define DECLARE_TABLES_COLL(_hf_type_, _hf_description_)                               \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##CollBases, "HF" _hf_description_ "COLLBASE", \
                           o2::soa::Index<>,                                           \
                           collision::PosX,                                            \
                           collision::PosY,                                            \
                           collision::PosZ,                                            \
                           collision::NumContrib,                                      \
                           hf_coll_base::CentFT0A,                                     \
                           hf_coll_base::CentFT0C,                                     \
                           hf_coll_base::CentFT0M,                                     \
                           hf_coll_base::CentFV0A,                                     \
                           hf_coll_base::MultZeqNTracksPV,                             \
                           o2::soa::Marker<Marker##_hf_type_>);                        \
                                                                                       \
  using Hf##_hf_type_##CollBase = Hf##_hf_type_##CollBases::iterator;                  \
  using StoredHf##_hf_type_##CollBase = StoredHf##_hf_type_##CollBases::iterator;      \
                                                                                       \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##CollIds, "HF" _hf_description_ "COLLID",     \
                           hf_cand::CollisionId,                                       \
                           o2::soa::Marker<Marker##_hf_type_>);

// Declares the base table with MC collisions (McCollBases) and joinable tables (McCollIds, McRCollIds).
#define DECLARE_TABLES_MCCOLL(_hf_type_, _hf_description_, _hf_namespace_)                                          \
  namespace hf_mc_coll                                                                                              \
  {                                                                                                                 \
  namespace der_##_hf_namespace_                                                                                    \
  {                                                                                                                 \
    DECLARE_SOA_ARRAY_INDEX_COLUMN_CUSTOM(Hf##_hf_type_##CollBase, hfCollBases, "HF" _hf_description_ "COLLBASES"); \
  }                                                                                                                 \
  }                                                                                                                 \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##McCollBases, "HF" _hf_description_ "MCCOLLBASE",                          \
                           o2::soa::Index<>,                                                                        \
                           mccollision::PosX,                                                                       \
                           mccollision::PosY,                                                                       \
                           mccollision::PosZ,                                                                       \
                           cent::CentFT0M,                                                                          \
                           o2::soa::Marker<Marker##_hf_type_>);                                                     \
                                                                                                                    \
  using Hf##_hf_type_##McCollBase = Hf##_hf_type_##McCollBases::iterator;                                           \
  using StoredHf##_hf_type_##McCollBase = StoredHf##_hf_type_##McCollBases::iterator;                               \
                                                                                                                    \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##McCollIds, "HF" _hf_description_ "MCCOLLID",                              \
                           hf_mc_coll::McCollisionId,                                                               \
                           o2::soa::Marker<Marker##_hf_type_>);                                                     \
                                                                                                                    \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##McRCollIds, "HF" _hf_description_ "MCRCOLLID",                            \
                           hf_mc_coll::der_##_hf_namespace_::Hf##_hf_type_##CollBaseIds);

// ================
// Candidate tables
// ================

namespace hf_cand_base
{
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

namespace hf_mc_particle
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);                 //! MC collision of this particle
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);                   //! MC particle
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! flag for generator level matching
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, generator level
} // namespace hf_mc_particle

// Declares the base table with candidates (Bases).
#define DECLARE_TABLE_CAND_BASE(_hf_type_, _hf_description_, _hf_namespace_)                                 \
  namespace hf_cand_base                                                                                     \
  {                                                                                                          \
  namespace der_##_hf_namespace_                                                                             \
  {                                                                                                          \
    DECLARE_SOA_INDEX_COLUMN_CUSTOM(Hf##_hf_type_##CollBase, hfCollBase, "HF" _hf_description_ "COLLBASES"); \
  }                                                                                                          \
  }                                                                                                          \
                                                                                                             \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##Bases, "HF" _hf_description_ "BASE",                               \
                           o2::soa::Index<>,                                                                 \
                           hf_cand_base::der_##_hf_namespace_::Hf##_hf_type_##CollBaseId,                    \
                           hf_cand_base::Pt,                                                                 \
                           hf_cand_base::Eta,                                                                \
                           hf_cand_base::Phi,                                                                \
                           hf_cand_base::M,                                                                  \
                           hf_cand_base::Y,                                                                  \
                           hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,                            \
                           hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,                            \
                           hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,                            \
                           hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,                             \
                           o2::soa::Marker<Marker##_hf_type_>);

// Declares the table with global indices for 2-prong candidates (Ids).
#define DECLARE_TABLE_CAND_ID_2P(_hf_type_, _hf_description_)              \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##Ids, "HF" _hf_description_ "ID", \
                           hf_cand::CollisionId,                           \
                           hf_track_index::Prong0Id,                       \
                           hf_track_index::Prong1Id,                       \
                           o2::soa::Marker<Marker##_hf_type_>);

// Declares the table with global indices for 3-prong candidates (Ids).
#define DECLARE_TABLE_CAND_ID_3P(_hf_type_, _hf_description_)              \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##Ids, "HF" _hf_description_ "ID", \
                           hf_cand::CollisionId,                           \
                           hf_track_index::Prong0Id,                       \
                           hf_track_index::Prong1Id,                       \
                           hf_track_index::Prong2Id,                       \
                           o2::soa::Marker<Marker##_hf_type_>);

// Declares the table with candidate selection flags (Sels).
#define DECLARE_TABLE_CAND_SEL(_hf_type_, _hf_description_)                  \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##Sels, "HF" _hf_description_ "SEL", \
                           hf_cand_sel::CandidateSelFlag,                    \
                           o2::soa::Marker<Marker##_hf_type_>);

// ================
// MC particle tables
// ================

// Declares the base table with MC particles (PBases).
#define DECLARE_TABLE_MCPARTICLE_BASE(_hf_type_, _hf_description_, _hf_namespace_)                                 \
  namespace hf_mc_particle                                                                                         \
  {                                                                                                                \
  namespace der_##_hf_namespace_                                                                                   \
  {                                                                                                                \
    DECLARE_SOA_INDEX_COLUMN_CUSTOM(Hf##_hf_type_##McCollBase, hfMcCollBase, "HF" _hf_description_ "MCCOLLBASES"); \
  }                                                                                                                \
  }                                                                                                                \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##PBases, "HF" _hf_description_ "PBASE",                                   \
                           o2::soa::Index<>,                                                                       \
                           hf_mc_particle::der_##_hf_namespace_::Hf##_hf_type_##McCollBaseId,                      \
                           hf_cand_base::Pt,                                                                       \
                           hf_cand_base::Eta,                                                                      \
                           hf_cand_base::Phi,                                                                      \
                           hf_cand_base::Y,                                                                        \
                           hf_mc_particle::FlagMcMatchGen,                                                         \
                           hf_mc_particle::OriginMcGen,                                                            \
                           hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,                                  \
                           hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,                                  \
                           hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,                                  \
                           hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,                                   \
                           o2::soa::Marker<Marker##_hf_type_>);

// Declares the table with global indices for MC particles (PIds).
#define DECLARE_TABLE_MCPARTICLE_ID(_hf_type_, _hf_description_)             \
  DECLARE_SOA_TABLE_STAGED(Hf##_hf_type_##PIds, "HF" _hf_description_ "PID", \
                           hf_mc_particle::McCollisionId,                    \
                           hf_mc_particle::McParticleId,                     \
                           o2::soa::Marker<Marker##_hf_type_>);

// ================
// Helper macros for combinations
// ================

#define DECLARE_TABLES_COMMON(_hf_type_, _hf_description_, _hf_namespace_)   \
  DECLARE_TABLES_COLL(_hf_type_, _hf_description_)                           \
  DECLARE_TABLES_MCCOLL(_hf_type_, _hf_description_, _hf_namespace_)         \
  DECLARE_TABLE_CAND_BASE(_hf_type_, _hf_description_, _hf_namespace_)       \
  DECLARE_TABLE_CAND_SEL(_hf_type_, _hf_description_)                        \
  DECLARE_TABLE_MCPARTICLE_BASE(_hf_type_, _hf_description_, _hf_namespace_) \
  DECLARE_TABLE_MCPARTICLE_ID(_hf_type_, _hf_description_)

#define DECLARE_TABLES_2P(_hf_type_, _hf_description_, _hf_namespace_, _marker_number_) \
  constexpr uint Marker##_hf_type_ = _marker_number_;                                   \
  DECLARE_TABLES_COMMON(_hf_type_, _hf_description_, _hf_namespace_)                    \
  DECLARE_TABLE_CAND_ID_2P(_hf_type_, _hf_description_)

#define DECLARE_TABLES_3P(_hf_type_, _hf_description_, _hf_namespace_, _marker_number_) \
  constexpr uint Marker##_hf_type_ = _marker_number_;                                   \
  DECLARE_TABLES_COMMON(_hf_type_, _hf_description_, _hf_namespace_)                    \
  DECLARE_TABLE_CAND_ID_3P(_hf_type_, _hf_description_)

// ================
// Declarations of common tables for individual species
// ================

DECLARE_TABLES_2P(D0, "D0", d0, 2);
DECLARE_TABLES_3P(Lc, "LC", lc, 3);
DECLARE_TABLES_3P(Dplus, "DP", dplus, 4);
DECLARE_TABLES_3P(Bplus, "BP", bplus, 5);

// ================
// Additional species-specific candidate tables
// ================

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

// ----------------
// D0
// ----------------

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

DECLARE_SOA_TABLE_STAGED(HfD0Mls, "HFD0ML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScores,
                         o2::soa::Marker<MarkerD0>);

DECLARE_SOA_TABLE_STAGED(HfD0Mcs, "HFD0MC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         o2::soa::Marker<MarkerD0>);

// ----------------
// B+
// ----------------

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

DECLARE_SOA_TABLE_STAGED(HfBplusMcs, "HFBPMC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         o2::soa::Marker<MarkerBplus>);

// ----------------
// Lc
// ----------------

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE_STAGED(HfLcPars, "HFLCPAR", //! Table with candidate properties used for selection
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
                         o2::soa::Marker<MarkerLc>);

DECLARE_SOA_TABLE_STAGED(HfLcParEs, "HFLCPARE", //! Table with additional candidate properties used for selection
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
                         o2::soa::Marker<MarkerLc>);

DECLARE_SOA_TABLE_STAGED(HfLcMls, "HFLCML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScores,
                         o2::soa::Marker<MarkerLc>);

DECLARE_SOA_TABLE_STAGED(HfLcMcs, "HFLCMC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         hf_cand_mc::IsCandidateSwapped,
                         o2::soa::Marker<MarkerLc>);

// ----------------
// D+
// ----------------

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE_STAGED(HfDplusPars, "HFDPPAR", //! Table with candidate properties used for selection
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
                         hf_cand_par::NSigTofPi0,
                         hf_cand_par::NSigTpcTofPi0,
                         hf_cand_par::NSigTpcKa1,
                         hf_cand_par::NSigTofKa1,
                         hf_cand_par::NSigTpcTofKa1,
                         hf_cand_par::NSigTpcPi2,
                         hf_cand_par::NSigTofPi2,
                         hf_cand_par::NSigTpcTofPi2,
                         o2::soa::Marker<MarkerDplus>);

DECLARE_SOA_TABLE_STAGED(HfDplusParEs, "HFDPPARE", //! Table with additional candidate properties used for selection
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
                         o2::soa::Marker<MarkerDplus>);

DECLARE_SOA_TABLE_STAGED(HfDplusMls, "HFDPML", //! Table with candidate selection ML scores
                         hf_cand_mc::MlScores,
                         o2::soa::Marker<MarkerDplus>);

DECLARE_SOA_TABLE_STAGED(HfDplusMcs, "HFDPMC", //! Table with MC candidate info
                         hf_cand_mc::FlagMcMatchRec,
                         hf_cand_mc::OriginMcRec,
                         hf_cand_mc::IsCandidateSwapped, // useless
                         hf_cand_mc::FlagMcDecayChanRec,
                         o2::soa::Marker<MarkerDplus>);

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLES_H_
