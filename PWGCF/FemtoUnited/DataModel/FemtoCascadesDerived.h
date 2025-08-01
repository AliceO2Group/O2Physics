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

/// \file FemtoCascadesDerived.h
/// \brief cascade tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_

#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"

#include "Framework/ASoA.h"

namespace o2::aod
{
namespace femtocascades
{
// columns for cascade bit masks
DECLARE_SOA_COLUMN(Mask, mask, femtodatatypes::CascadeMaskType); //! Bitmask for cascade selections

// columns for cascad debug information
DECLARE_SOA_COLUMN(MassXi, massXi, float);                         //! Mass of xi
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);                   //! Mass of omega
DECLARE_SOA_COLUMN(CascadeCosPa, cascadeCosPa, float);             //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(CascadeDauDca, cascadeDauDca, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(CascadeTransRadius, cascadeTransRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaCosPa, lambdaCosPa, float);               //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(LambdaDauDca, lambdaDauDca, float);             //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(LambdaTransRadius, lambdaTransRadius, float);   //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaDcaToPv, lambdaDcaToPv, float);           //! Lambda transvers radius

// id columns for bachelor
// following same style as strangeness tables were we do not store the id of the lambda, but its daughters
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int32_t, FUTracks, "_Bachelor"); //! bachelor id

} // namespace femtocascades

// table for basic vzero information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUXis_001, "FUXIS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtocascades::BachelorId,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUXis = FUXis_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUXiMasks_001, "FUXIMASKS", 1,
                                   femtocascades::Mask);
using FUXiMasks = FUXiMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUOmegas_001, "FUOMEGAS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::SignedPt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtocascades::BachelorId,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUOmegas = FUOmegas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUOmegaMasks_001, "FUOMEGAMASKS", 1,
                                   femtocascades::Mask);
using FUOmegaMasks = FUOmegaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUXiExtras_001, "FUXIEXTRAS", 1,
                                   femtocascades::MassOmega,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv);
using FUXiExtras = FUXiExtras_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUOmegaExtras_001, "FUOMEGAEXTRAS", 1,
                                   femtocascades::MassXi,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv);
using FUOmegaExtras = FUOmegaExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
