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

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace o2::aod
{
// Basic collision properties
namespace hf_coll_base
{
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int8_t); //! collision rejection flag
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);            //! FT0M multiplicity
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);            //! FT0A centrality percentile
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);            //! FT0C centrality percentile
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);            //! FT0M centrality percentile
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);            //! FT0A centrality percentile
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);            //! FDDM centrality percentile
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float); //! z-equalised barrel multiplicity
} // namespace hf_coll_base

// D0

DECLARE_SOA_TABLE(HfD0CollBases, "AOD", "HFD0COLLBASE", //! Table with basic collision info
                  o2::soa::Index<>,
                  collision::NumContrib,
                  hf_coll_base::IsEventReject,
                  bc::RunNumber);

using HfD0CollBase = HfD0CollBases::iterator;

DECLARE_SOA_TABLE(StoredHfD0CollBases, "AOD1", "HFD0COLLBASE", //! Table with basic collision info (stored version)
                  o2::soa::Index<>,
                  collision::NumContrib,
                  hf_coll_base::IsEventReject,
                  bc::RunNumber,
                  soa::Marker<1>);

using StoredHfD0CollBase = StoredHfD0CollBases::iterator;

DECLARE_SOA_TABLE(HfD0CollIds, "AOD", "HFD0COLLID", //! Table with global indices for collisions
                  hf_cand::CollisionId);

DECLARE_SOA_TABLE(StoredHfD0CollIds, "AOD1", "HFD0COLLID", //! Table with global indices for collisions (stored version)
                  hf_cand::CollisionId,
                  soa::Marker<1>);

// 3-prong decays

DECLARE_SOA_TABLE(Hf3PCollBases, "AOD", "HF3PCOLLBASE", //! Table with basic collision info
                  o2::soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::NumContrib,
                  hf_coll_base::CentFT0A,
                  hf_coll_base::CentFT0C,
                  hf_coll_base::CentFT0M,
                  hf_coll_base::CentFV0A,
                  hf_coll_base::MultZeqNTracksPV);
// hf_coll_base::IsEventReject,
// bc::RunNumber,

using Hf3PCollBase = Hf3PCollBases::iterator;

DECLARE_SOA_TABLE(StoredHf3PCollBases, "AOD1", "HF3PCOLLBASE", //! Table with basic collision info (stored version)
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
                  soa::Marker<1>);

using StoredHf3PCollBase = StoredHf3PCollBases::iterator;

DECLARE_SOA_TABLE(Hf3PCollIds, "AOD", "HF3PCOLLID", //! Table with global indices for collisions
                  hf_cand::CollisionId,
                  soa::Marker<2>);

DECLARE_SOA_TABLE(StoredHf3PCollIds, "AOD1", "HF3PCOLLID", //! Table with global indices for collisions (stored version)
                  hf_cand::CollisionId,
                  soa::Marker<3>);

// Basic candidate properties
namespace hf_cand_base
{
DECLARE_SOA_INDEX_COLUMN(HfD0CollBase, hfD0CollBase);             //! collision index pointing to the derived collision table for D0 candidates
DECLARE_SOA_INDEX_COLUMN(StoredHfD0CollBase, storedHfD0CollBase); //! collision index pointing to the derived collision table for D0 candidates
DECLARE_SOA_INDEX_COLUMN(Hf3PCollBase, hf3PCollBase);             //! collision index pointing to the derived collision table for 3-prong candidates
DECLARE_SOA_INDEX_COLUMN(StoredHf3PCollBase, storedHf3PCollBase); //! collision index pointing to the derived collision table for 3-prong candidates
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);               //! MC collision of this particle
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);                 //! MC particle
DECLARE_SOA_COLUMN(Eta, eta, float);                              //! pseudorapidity
DECLARE_SOA_COLUMN(M, m, float);                                  //! invariant mass
DECLARE_SOA_COLUMN(Phi, phi, float);                              //! azimuth
DECLARE_SOA_COLUMN(Pt, pt, float);                                //! transverse momentum

namespace functions_pt_eta_phi
{
/// px as a function of pT, phi
/// \todo Move to RecoDecay
template <typename TPt, typename TPhi>
auto px(TPt pt, TPhi phi)
{
  return pt * std::cos(phi);
}

/// py as a function of pT, phi
/// \todo Move to RecoDecay
template <typename TPt, typename TPhi>
auto py(TPt pt, TPhi phi)
{
  return pt * std::sin(phi);
}

/// pz as a function of pT, eta
/// \todo Move to RecoDecay
template <typename TPt, typename TEta>
auto pz(TPt pt, TEta eta)
{
  return pt * std::sinh(eta);
}

/// p as a function of pT, eta
/// \todo Move to RecoDecay
template <typename TPt, typename TEta>
auto p(TPt pt, TEta eta)
{
  return pt * std::cosh(eta);
}

/// Rapidity as a function of pT, eta, mass
/// \todo Move to RecoDecay
template <typename TPt, typename TEta, typename TM>
auto y(TPt pt, TEta eta, TM m)
{
  return std::log((RecoDecay::sqrtSumOfSquares(m, pt * std::cosh(eta)) + pt * std::sinh(eta)) / RecoDecay::sqrtSumOfSquares(m, pt));
}

/// Energy as a function of pT, eta, mass
/// \todo Move to RecoDecay
template <typename TPt, typename TEta, typename TM>
auto e(TPt pt, TEta eta, TM m)
{
  return RecoDecay::sqrtSumOfSquares(m, p(pt, eta));
}
} // namespace functions_pt_eta_phi

namespace d0
{
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! D0 rapidity
                           [](float pt, float eta) -> float { return functions_pt_eta_phi::y(pt, eta, o2::constants::physics::MassD0); });
}
namespace lc
{
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Lambda_c rapidity
                           [](float pt, float eta) -> float { return functions_pt_eta_phi::y(pt, eta, o2::constants::physics::MassLambdaCPlus); });
}
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! px
                           [](float pt, float phi) -> float { return functions_pt_eta_phi::px(pt, phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! py
                           [](float pt, float phi) -> float { return functions_pt_eta_phi::py(pt, phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! px
                           [](float pt, float eta) -> float { return functions_pt_eta_phi::pz(pt, eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! momentum
                           [](float pt, float eta) -> float { return functions_pt_eta_phi::p(pt, eta); });
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

// Candidate selection flags
namespace hf_cand_sel
{
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t); //! bitmap of the selected candidate type
}

// MC flags
namespace hf_cand_mc
{
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);         //! flag for reconstruction level matching
DECLARE_SOA_COLUMN(FlagMcMatchGen, flagMcMatchGen, int8_t);         //! flag for generator level matching
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);               //! particle origin, reconstruction level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);               //! particle origin, generator level
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t); //! swapping of the prongs order
DECLARE_SOA_COLUMN(FlagMcDecayChanRec, flagMcDecayChanRec, int8_t); //! resonant decay channel flag, reconstruction level
DECLARE_SOA_COLUMN(FlagMcDecayChanGen, flagMcDecayChanGen, int8_t); //! resonant decay channel flag, generator level
DECLARE_SOA_COLUMN(MlScoreBkg, mlScoreBkg, float);                  //! ML score for background class
DECLARE_SOA_COLUMN(MlScorePrompt, mlScorePrompt, float);            //! ML score for prompt class
DECLARE_SOA_COLUMN(MlScoreNonPrompt, mlScoreNonPrompt, float);      //! ML score for non-prompt class
} // namespace hf_cand_mc

// D0

DECLARE_SOA_TABLE(HfD0Bases, "AOD", "HFD0BASE", //! Table with basic candidate properties used in the analyses
                  o2::soa::Index<>,
                  hf_cand_base::HfD0CollBaseId,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::d0::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>);

DECLARE_SOA_TABLE(StoredHfD0Bases, "AOD1", "HFD0BASE", //! Table with basic candidate properties used in the analyses (stored version)
                  o2::soa::Index<>,
                  hf_cand_base::HfD0CollBaseId,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::d0::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                  soa::Marker<1>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(HfD0Pars, "AOD", "HFD0PAR", //! Table with candidate properties used for selection
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

DECLARE_SOA_TABLE(StoredHfD0Pars, "AOD1", "HFD0PAR", //! Table with candidate properties used for selection (stored version)
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
                  hf_cand_par::ImpactParameterProduct,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0ParEs, "AOD", "HFD0PARE", //! Table with additional candidate properties used for selection
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

DECLARE_SOA_TABLE(StoredHfD0ParEs, "AOD1", "HFD0PARE", //! Table with additional candidate properties used for selection (stored version)
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
                  hf_cand_par::Ct,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0Sels, "AOD", "HFD0SEL", //! Table with candidate selection flags
                  hf_cand_sel::CandidateSelFlag);

DECLARE_SOA_TABLE(StoredHfD0Sels, "AOD1", "HFD0SEL", //! Table with candidate selection flags (stored version)
                  hf_cand_sel::CandidateSelFlag,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0Ids, "AOD", "HFD0ID", //! Table with global indices for candidates
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);

DECLARE_SOA_TABLE(StoredHfD0Ids, "AOD1", "HFD0ID", //! Table with global indices for candidates (stored version)
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0Mcs, "AOD", "HFD0MC", //! Table with MC candidate info
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec);

DECLARE_SOA_TABLE(StoredHfD0Mcs, "AOD1", "HFD0MC", //! Table with MC candidate info (stored version)
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0PBases, "AOD", "HFD0PBASE", //! Table with MC particle info
                  o2::soa::Index<>,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_mc::FlagMcMatchGen,
                  hf_cand_mc::OriginMcGen,
                  hf_cand_base::d0::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>);

DECLARE_SOA_TABLE(StoredHfD0PBases, "AOD1", "HFD0PBASE", //! Table with MC particle info (stored version)
                  o2::soa::Index<>,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_mc::FlagMcMatchGen,
                  hf_cand_mc::OriginMcGen,
                  hf_cand_base::d0::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(HfD0PIds, "AOD", "HFD0PID", //! Table with global indices for MC particles
                  hf_cand_base::McCollisionId,
                  hf_cand_base::McParticleId);

DECLARE_SOA_TABLE(StoredHfD0PIds, "AOD1", "HFD0PID", //! Table with global indices for MC particles (stored version)
                  hf_cand_base::McCollisionId,
                  hf_cand_base::McParticleId,
                  soa::Marker<1>);

// 3-prong decays

DECLARE_SOA_TABLE(Hf3PBases, "AOD", "HF3PBASE", //! Table with basic candidate properties used in the analyses
                  o2::soa::Index<>,
                  hf_cand_base::Hf3PCollBaseId,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::lc::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>);

DECLARE_SOA_TABLE(StoredHf3PBases, "AOD1", "HF3PBASE", //! Table with basic candidate properties used in the analyses (stored version)
                  o2::soa::Index<>,
                  hf_cand_base::Hf3PCollBaseId,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_base::M,
                  hf_cand_base::lc::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                  soa::Marker<1>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(Hf3PPars, "AOD", "HF3PPAR", //! Table with candidate properties used for selection
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
                  hf_cand_par::NSigTpcTofPr2);

DECLARE_SOA_TABLE(StoredHf3PPars, "AOD1", "HF3PPAR", //! Table with candidate properties used for selection (stored version)
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
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_par::ImpactParameterNormalised0,
                  hf_cand_par::ImpactParameterNormalised1,
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
                  soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PParEs, "AOD", "HF3PPARE", //! Table with additional candidate properties used for selection
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
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
                  hf_cand_par::Ct);

DECLARE_SOA_TABLE(StoredHf3PParEs, "AOD1", "HF3PPARE", //! Table with additional candidate properties used for selection (stored version)
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength,
                  hf_cand::ErrorDecayLengthXY,
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
                  hf_cand_par::Ct,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PSels, "AOD", "HF3PSEL", //! Table with candidate selection flags
                  hf_cand_sel::CandidateSelFlag,
                  soa::Marker<2>);

DECLARE_SOA_TABLE(StoredHf3PSels, "AOD1", "HF3PSEL", //! Table with candidate selection flags (stored version)
                  hf_cand_sel::CandidateSelFlag,
                  soa::Marker<3>);

DECLARE_SOA_TABLE(Hf3PIds, "AOD", "HF3PID", //! Table with global indices for candidates
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_track_index::Prong2Id);

DECLARE_SOA_TABLE(StoredHf3PIds, "AOD1", "HF3PID", //! Table with global indices for candidates (stored version)
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_track_index::Prong2Id,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PMcs, "AOD", "HF3PMC", //! Table with MC candidate info
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec,
                  hf_cand_mc::IsCandidateSwapped);

DECLARE_SOA_TABLE(StoredHf3PMcs, "AOD1", "HF3PMC", //! Table with MC candidate info (stored version)
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec,
                  hf_cand_mc::IsCandidateSwapped,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PPBases, "AOD", "HF3PPBASE", //! Table with MC particle info
                  o2::soa::Index<>,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_mc::FlagMcMatchGen,
                  hf_cand_mc::OriginMcGen,
                  hf_cand_base::lc::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>);

DECLARE_SOA_TABLE(StoredHf3PPBases, "AOD1", "HF3PPBASE", //! Table with MC particle info (stored version)
                  o2::soa::Index<>,
                  hf_cand_base::Pt,
                  hf_cand_base::Eta,
                  hf_cand_base::Phi,
                  hf_cand_mc::FlagMcMatchGen,
                  hf_cand_mc::OriginMcGen,
                  hf_cand_base::lc::Y<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::Px<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Py<hf_cand_base::Pt, hf_cand_base::Phi>,
                  hf_cand_base::Pz<hf_cand_base::Pt, hf_cand_base::Eta>,
                  hf_cand_base::P<hf_cand_base::Pt, hf_cand_base::Eta>,
                  soa::Marker<1>);

DECLARE_SOA_TABLE(Hf3PPIds, "AOD", "HF3PPID", //! Table with global indices for MC particles
                  hf_cand_base::McCollisionId,
                  hf_cand_base::McParticleId,
                  soa::Marker<2>);

DECLARE_SOA_TABLE(StoredHf3PPIds, "AOD1", "HF3PPID", //! Table with global indices for MC particles (stored version)
                  hf_cand_base::McCollisionId,
                  hf_cand_base::McParticleId,
                  soa::Marker<3>);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLES_H_
