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
/// \file lambdaJetPolarizationIons.h
/// \brief Derived Data table for Jet-induced polarization analysis (HI)
/// \author Cicero Domenico Muncinelli (cicero.domenico.muncinelli@cern.ch)
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    cicero.domenico.muncinelli@cern.ch
//

#ifndef PWGLF_DATAMODEL_LAMBDAJETPOL_H_
#define PWGLF_DATAMODEL_LAMBDAJETPOL_H_

#include <Framework/ASoA.h>

namespace o2::aod
{

namespace lambdajetpol
{

// DECLARE_SOA_COLUMN(CollIdx, collIdx, uint64_t); // Using a regular SOA column instead of an index column for convenience
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Centrality, centrality, float);

DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, uint64_t);

DECLARE_SOA_COLUMN(V0Pt, v0Pt, float);
DECLARE_SOA_COLUMN(V0Eta, v0Eta, float);
DECLARE_SOA_COLUMN(V0Phi, v0Phi, float);

DECLARE_SOA_COLUMN(IsLambda, isLambda, bool);
DECLARE_SOA_COLUMN(IsAntiLambda, isAntiLambda, bool);
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);
DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float);

DECLARE_SOA_COLUMN(PosPt, posPt, float);
DECLARE_SOA_COLUMN(PosEta, posEta, float);
DECLARE_SOA_COLUMN(PosPhi, posPhi, float);
DECLARE_SOA_COLUMN(NegPt, negPt, float);
DECLARE_SOA_COLUMN(NegEta, negEta, float);
DECLARE_SOA_COLUMN(NegPhi, negPhi, float);

// (TODO: add dynamic columns with jet px, py, pz)
// (TODO: add leading particle as one of the columns!)

} // namespace lambdajetpol

DECLARE_SOA_TABLE(RingJets, "AOD", "RINGJETS", // Renamed to follow convention on "s" at the end of table name.
                  lambdajetpol::CollisionId, // Changed to an internal O2 index, slightly different from usual o2::soa::Index<> though
                  lambdajetpol::JetPt,
                  lambdajetpol::JetEta,
                  lambdajetpol::JetPhi,
                  lambdajetpol::JetNConstituents);

DECLARE_SOA_TABLE(RingLaV0s, "AOD", "RINGLAV0S", // Had to write this in a shorter form because the derived data did not accept long names
                  lambdajetpol::CollisionId,
                  lambdajetpol::V0Pt,
                  lambdajetpol::V0Eta,
                  lambdajetpol::V0Phi,
                  lambdajetpol::IsLambda,
                  lambdajetpol::IsAntiLambda,
                  lambdajetpol::MassLambda,
                  lambdajetpol::MassAntiLambda,
                  lambdajetpol::PosPt,
                  lambdajetpol::PosEta,
                  lambdajetpol::PosPhi,
                  lambdajetpol::NegPt,
                  lambdajetpol::NegEta,
                  lambdajetpol::NegPhi);

DECLARE_SOA_TABLE(RingCollisions, "AOD", "RINGCOLLISIONS",
                  lambdajetpol::CollisionId,
                  lambdajetpol::Centrality);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LAMBDAJETPOL_H_