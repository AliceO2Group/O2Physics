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
namespace hf_cand_analysis
{
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
}
namespace hf_cand_param
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
}
namespace hf_cand_flag
{
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
}
namespace hf_index
{
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate2P, candidate2p, int, HfCand2Prong, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate3P, candidate3p, int, HfCand3Prong, "_0");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
}
namespace hf_collision
{
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int8_t);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(MultZeqFT0A, multZeqFT0A, float);
DECLARE_SOA_COLUMN(MultZeqFT0C, multZeqFT0C, float);
DECLARE_SOA_COLUMN(MultZeqFV0A, multZeqFV0A, float);
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
}
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLES_H_
