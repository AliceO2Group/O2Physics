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

///
/// \file   OTFMCParticle.h
/// \author Jesper Karlsson Gumprecht
/// \since  16/12/2025
/// \brief  Redefinition of the mcparticles table specifically for the fast sim
///

#ifndef ALICE3_DATAMODEL_OTFMCPARTICLE_H_
#define ALICE3_DATAMODEL_OTFMCPARTICLE_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

namespace otfmcparticle
{
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(IsAlive, isAlive, bool);
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, bool);

} // namespace otfmcparticle

DECLARE_SOA_TABLE_FULL(McParticlesWithDau, "McParticlesWithDau", "AOD", "MCPARTICLEWITHDAU",
                       o2::soa::Index<>,
                       mcparticle::McCollisionId,
                       mcparticle::PdgCode,
                       mcparticle::StatusCode,
                       mcparticle::Flags,
                       mcparticle::MothersIds,
                       mcparticle::DaughtersIdSlice,
                       mcparticle::Weight,
                       mcparticle::Px,
                       mcparticle::Py,
                       mcparticle::Pz,
                       mcparticle::E,
                       mcparticle::Vx,
                       mcparticle::Vy,
                       mcparticle::Vz,
                       mcparticle::Vt,
                       otfmcparticle::Phi,
                       otfmcparticle::Eta,
                       otfmcparticle::Pt,
                       otfmcparticle::P,
                       otfmcparticle::Y,
                       otfmcparticle::IsAlive,
                       otfmcparticle::IsPrimary,
                       mcparticle::PVector<mcparticle::Px, mcparticle::Py, mcparticle::Pz>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetHepMCStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using McParticleWithDau = McParticlesWithDau::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMCPARTICLE_H_
