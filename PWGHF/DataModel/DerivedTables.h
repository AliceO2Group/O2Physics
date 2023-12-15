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

/// \file DerivedTables.h.h
/// \brief Definitions of derived tables produced by derived-data creators

#ifndef PWGHF_DATAMODEL_DERIVEDTABLES_H_
#define PWGHF_DATAMODEL_DERIVEDTABLES_H_

#include "Framework/ASoA.h"

namespace o2::aod
{
namespace hf_cand_index
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
} // namespace hf_cand_index

namespace hf_coll_base
{
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int8_t);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(MultZeqFT0A, multZeqFT0A, float);
DECLARE_SOA_COLUMN(MultZeqFT0C, multZeqFT0C, float);
DECLARE_SOA_COLUMN(MultZeqFV0A, multZeqFV0A, float);
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace hf_coll_base

// Table with basic collision info
DECLARE_SOA_TABLE(HfD0CollBases, "AOD", "HFD0COLLBASE",
                  o2::soa::Index<>,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_coll_base::IsEventReject,
                  hf_coll_base::RunNumber);

using HfD0CollBase = HfD0CollBases::iterator;

// Table with global indices for collisions
DECLARE_SOA_TABLE(HfD0CollIds, "AOD", "HFD0COLLID",
                  hf_cand_index::CollisionId);

namespace hf_cand_base
{
DECLARE_SOA_INDEX_COLUMN(HfD0CollBase, hfD0CollBase);  //! Collision index
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
} // namespace hf_cand_base

namespace hf_cand_par
{
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(PProng2, pProng2, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
// TOF
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);
DECLARE_SOA_COLUMN(NSigTofPr1, nSigTofPr1, float);
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);
// TPC
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);
DECLARE_SOA_COLUMN(NSigTpcPr1, nSigTpcPr1, float);
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);
// TPC+TOF
DECLARE_SOA_COLUMN(NSigTpcTofKa0, nSigTpcTofKa0, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa2, nSigTpcTofKa2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr0, nSigTpcTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr1, nSigTpcTofPr1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr2, nSigTpcTofPr2, float);
} // namespace hf_cand_par

namespace hf_cand_sel
{
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
}

namespace hf_cand_mc
{
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(MlScoreBkg, mlScoreBkg, float);             //! ML score for background class
DECLARE_SOA_COLUMN(MlScorePrompt, mlScorePrompt, float);       //! ML score for prompt class
DECLARE_SOA_COLUMN(MlScoreNonPrompt, mlScoreNonPrompt, float); //! ML score for non-prompt class
} // namespace hf_cand_mc

/// Table with basic candidate properties used in the analyses
// candidates for removal:
// E
DECLARE_SOA_TABLE(HfD0Bases, "AOD", "HFD0BASE",
                  o2::soa::Index<>,
                  hf_cand_base::HfD0CollBaseId,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::P,
                  hf_cand_base::E,
                  hf_cand_base::Y);

/// Table with candidate properties used for selection
// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep P, Pt, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(HfD0Pars, "AOD", "HFD0PAR",
                  hf_cand::Chi2PCA,
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
                  hf_cand_par::NSigTpcPi0,
                  hf_cand_par::NSigTpcKa0,
                  hf_cand_par::NSigTofPi0,
                  hf_cand_par::NSigTofKa0,
                  hf_cand_par::NSigTpcTofPi0,
                  hf_cand_par::NSigTpcTofKa0,
                  hf_cand_par::NSigTpcPi1,
                  hf_cand_par::NSigTpcKa1,
                  hf_cand_par::NSigTofPi1,
                  hf_cand_par::NSigTofKa1,
                  hf_cand_par::NSigTpcTofPi1,
                  hf_cand_par::NSigTpcTofKa1,
                  hf_cand_par::Cpa,
                  hf_cand_par::CpaXY,
                  hf_cand_par::MaxNormalisedDeltaIP,
                  hf_cand_par::ImpactParameterProduct);

/// Table with additional candidate properties used for selection
DECLARE_SOA_TABLE(HfD0ParEs, "AOD", "HFD0PARE",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
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
                  hf_cand_par::Ct);

// Table with candidate selection flags
DECLARE_SOA_TABLE(HfD0Sels, "AOD", "HFD0SEL",
                  hf_cand_sel::CandidateSelFlag);

// Table with global indices for candidates
DECLARE_SOA_TABLE(HfD0Ids, "AOD", "HFD0ID",
                  hf_track_index::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);

// Table with candidate MC info
DECLARE_SOA_TABLE(HfD0Mcs, "AOD", "HFD0MC",
                  hf_cand_mc::FlagMc,
                  hf_cand_mc::OriginMcRec);

// Table with MC particle info
DECLARE_SOA_TABLE(HfD0Ps, "AOD", "HFD0P",
                  o2::soa::Index<>,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::Y,
                  hf_cand_mc::FlagMc,
                  hf_cand_mc::OriginMcGen);

// Table with global indices for MC particles
DECLARE_SOA_TABLE(HfD0PIds, "AOD", "HFD0PID",
                  hf_cand_index::McCollisionId,
                  hf_cand_index::McParticleId);

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLES_H_
