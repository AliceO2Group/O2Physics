// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoV0sDerived.h
/// \brief v0 tables tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_

#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

namespace o2::aod
{
namespace femtov0s
{
// columns for bit masks
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::V0MaskType); //! Bitmask for v0 selections

// columns for debug information
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);         //! Mass of Lambda
DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float); //! Mass of AntiLambda
DECLARE_SOA_COLUMN(MassK0short, massK0short, float);       //! Mass of K0short
DECLARE_SOA_COLUMN(CosPa, cosPa, float);                   //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(DauDca, dauDca, float);                 //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);       //! Lambda transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);           //! x coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);           //! y coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);           //! z coordinate of Lambda decay vertex
DECLARE_SOA_DYNAMIC_COLUMN(DecayVtx, decayVtx,             //! distance of decay vertex from nominal interaction point
                           [](float vtxX, float vtxY, float vtxZ) -> float {
                             return std::hypot(vtxX, vtxY, vtxZ);
                           });

// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FUTracks, "_PosDau"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FUTracks, "_NegDau"); //!

} // namespace femtov0s

// table for basic lambda information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdas_001, "FULAMBDAS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId, // use sign to differentiate between lambda (+1) and antilambda (-1)
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass, // mass of the lambda/antilambda depending on the sign of the pt
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FULambdas = FULambdas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaMasks_001, "FULAMBDAMASKS", 1,
                                   femtov0s::Mask);
using FULambdaMasks = FULambdaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaExtras_001, "FULAMBDAEXTRAS", 1,
                                   femtobase::stored::MassAnti, // put mass of antiparticle, i.e. antilambda mass for lambdas and vice versa
                                   femtov0s::MassK0short,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
                                   femtov0s::TransRadius,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::DecayVtx<femtov0s::DecayVtxX, femtov0s::DecayVtxY, femtov0s::DecayVtxZ>);

using FULambdaExtras = FULambdaExtras_001;

// table for basic k0short information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUK0shorts_001, "FUK0SHORTS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUK0shorts = FUK0shorts_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUK0shortMasks_001, "FUK0SHORTMASKS", 1,
                                   femtov0s::Mask);
using FUK0shortMasks = FUK0shortMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUK0shortExtras_001, "FUK0SHORTEXTRAS", 1,
                                   femtov0s::MassLambda,
                                   femtov0s::MassAntiLambda,
                                   femtov0s::CosPa,
                                   femtov0s::DauDca,
                                   femtov0s::TransRadius,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::DecayVtx<femtov0s::DecayVtxX, femtov0s::DecayVtxY, femtov0s::DecayVtxZ>);

using FUK0shortExtras = FUK0shortExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_
