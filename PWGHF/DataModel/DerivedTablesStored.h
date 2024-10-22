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

/// \file DerivedTablesStored.h
/// \brief Definitions of stored versions of derived tables produced by derived-data creators (defined in DerivedTables.h)
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \note Things to check when comparing with DerivedTables.h:
///       - Prefix "Stored" in table definitions
///       - Prefix "Stored" in names of index columns pointing to derived tables
///       - Suffix "Stored" in Marker name
///       - Prefix "der_stored_" in namespace names (if needed to avoid redefinitions in "der_")
///       - Origin AOD1

#ifndef PWGHF_DATAMODEL_DERIVEDTABLESSTORED_H_
#define PWGHF_DATAMODEL_DERIVEDTABLESSTORED_H_

#include "PWGHF/DataModel/DerivedTables.h"

namespace o2::aod
{
constexpr uint MarkerD0Stored = 10;
constexpr uint Marker3PStored = 20;

// ================
// Collision tables
// ================

// D0

DECLARE_SOA_TABLE(StoredHfD0CollBases, "AOD1", "HFD0COLLBASE", //! Table with basic collision info
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
                  soa::Marker<MarkerD0Stored>);

using StoredHfD0CollBase = StoredHfD0CollBases::iterator;

DECLARE_SOA_TABLE(StoredHfD0CollIds, "AOD1", "HFD0COLLID", //! Table with original global indices of collisions
                  hf_cand::CollisionId,
                  soa::Marker<MarkerD0Stored>);

// 3-prong decays

DECLARE_SOA_TABLE(StoredHf3PCollBases, "AOD1", "HF3PCOLLBASE", //! Table with basic collision info
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
                  soa::Marker<Marker3PStored>);

using StoredHf3PCollBase = StoredHf3PCollBases::iterator;

DECLARE_SOA_TABLE(StoredHf3PCollIds, "AOD1", "HF3PCOLLID", //! Table with original global indices of collisions
                  hf_cand::CollisionId,
                  soa::Marker<Marker3PStored>);

// ===================
// MC collision tables
// ===================

// DO

DECLARE_SOA_TABLE(StoredHfD0McCollBases, "AOD1", "HFD0MCCOLLBASE", //! Table with basic MC collision info
                  o2::soa::Index<>,
                  mccollision::PosX,
                  mccollision::PosY,
                  mccollision::PosZ,
                  soa::Marker<MarkerD0Stored>);

using StoredHfD0McCollBase = StoredHfD0McCollBases::iterator;

DECLARE_SOA_TABLE(StoredHfD0McCollIds, "AOD1", "HFD0MCCOLLID", //! Table with original global indices of MC collisions
                  hf_mc_coll::McCollisionId,
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0McRCollIds, "AOD1", "HFD0MCRCOLLID", //! Table with indices pointing to the derived reconstructed-collision table
                  hf_mc_coll::der_d0::HfD0CollBaseIds,
                  soa::Marker<MarkerD0Stored>);

// 3-prong decays

DECLARE_SOA_TABLE(StoredHf3PMcCollBases, "AOD1", "HF3PMCCOLLBASE", //! Table with basic MC collision info
                  o2::soa::Index<>,
                  mccollision::PosX,
                  mccollision::PosY,
                  mccollision::PosZ,
                  soa::Marker<Marker3PStored>);

using StoredHf3PMcCollBase = StoredHf3PMcCollBases::iterator;

DECLARE_SOA_TABLE(StoredHf3PMcCollIds, "AOD1", "HF3PMCCOLLID", //! Table with original global indices of MC collisions
                  hf_mc_coll::McCollisionId,
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PMcRCollIds, "AOD1", "HF3PMCRCOLLID", //! Table with indices pointing to the derived reconstructed-collision table
                  hf_mc_coll::der_3p::Hf3PCollBaseIds,
                  soa::Marker<Marker3PStored>);

// ================
// Candidate tables
// ================

// Basic candidate properties

// D0

DECLARE_SOA_TABLE(StoredHfD0Bases, "AOD1", "HFD0BASE", //! Table with basic candidate properties used in the analyses
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
                  soa::Marker<MarkerD0Stored>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(StoredHfD0Pars, "AOD1", "HFD0PAR", //! Table with candidate properties used for selection
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
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0ParEs, "AOD1", "HFD0PARE", //! Table with additional candidate properties used for selection
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
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0Sels, "AOD1", "HFD0SEL", //! Table with candidate selection flags
                  hf_cand_sel::CandidateSelFlag,
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0Mls, "AOD1", "HFD0ML", //! Table with candidate selection ML scores
                  hf_cand_mc::MlScores,
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0Ids, "AOD1", "HFD0ID", //! Table with original global indices for candidates
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0Mcs, "AOD1", "HFD0MC", //! Table with MC candidate info
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec,
                  soa::Marker<MarkerD0Stored>);

// 3-prong decays

DECLARE_SOA_TABLE(StoredHf3PBases, "AOD1", "HF3PBASE", //! Table with basic candidate properties used in the analyses
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
                  soa::Marker<Marker3PStored>);

// candidates for removal:
// PxProng0, PyProng0, PzProng0,... (same for 1, 2), we can keep Pt, Eta, Phi instead
// XY: CpaXY, DecayLengthXY, ErrorDecayLengthXY
// normalised: DecayLengthNormalised, DecayLengthXYNormalised, ImpactParameterNormalised0
DECLARE_SOA_TABLE(StoredHf3PPars, "AOD1", "HF3PPAR", //! Table with candidate properties used for selection
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
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PParEs, "AOD1", "HF3PPARE", //! Table with additional candidate properties used for selection
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
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PSels, "AOD1", "HF3PSEL", //! Table with candidate selection flags
                  hf_cand_sel::CandidateSelFlag,
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PMls, "AOD1", "HF3PML", //! Table with candidate selection ML scores
                  hf_cand_mc::MlScores,
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PIds, "AOD1", "HF3PID", //! Table with original global indices for candidates
                  hf_cand::CollisionId,
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id,
                  hf_track_index::Prong2Id,
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PMcs, "AOD1", "HF3PMC", //! Table with MC candidate info
                  hf_cand_mc::FlagMcMatchRec,
                  hf_cand_mc::OriginMcRec,
                  hf_cand_mc::IsCandidateSwapped,
                  soa::Marker<Marker3PStored>);

// ==================
// MC particle tables
// ==================

// D0

DECLARE_SOA_TABLE(StoredHfD0PBases, "AOD1", "HFD0PBASE", //! Table with MC particle info
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
                  soa::Marker<MarkerD0Stored>);

DECLARE_SOA_TABLE(StoredHfD0PIds, "AOD1", "HFD0PID", //! Table with original global indices for MC particles
                  hf_mc_particle::McCollisionId,
                  hf_mc_particle::McParticleId,
                  soa::Marker<MarkerD0Stored>);

// 3-prong decays

DECLARE_SOA_TABLE(StoredHf3PPBases, "AOD1", "HF3PPBASE", //! Table with MC particle info
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
                  soa::Marker<Marker3PStored>);

DECLARE_SOA_TABLE(StoredHf3PPIds, "AOD1", "HF3PPID", //! Table with original global indices for MC particles
                  hf_mc_particle::McCollisionId,
                  hf_mc_particle::McParticleId,
                  soa::Marker<Marker3PStored>);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVEDTABLESSTORED_H_
