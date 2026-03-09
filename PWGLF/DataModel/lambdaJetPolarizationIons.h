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

#ifndef PWGLF_DATAMODEL_LAMBDAJETPOLARIZATIONIONS_H_
#define PWGLF_DATAMODEL_LAMBDAJETPOLARIZATIONIONS_H_

#include <Framework/ASoA.h>
#include <cmath>

namespace o2::aod
{

namespace lambdajetpol
{
// Collision information:
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);

// Jet (and jet proxies) information:
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, int);

DECLARE_SOA_COLUMN(LeadParticlePt, leadParticlePt, float);
DECLARE_SOA_COLUMN(LeadParticleEta, leadParticleEta, float);
DECLARE_SOA_COLUMN(LeadParticlePhi, leadParticlePhi, float);

// V0 information:
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

DECLARE_SOA_COLUMN(PosTPCNSigmaPr, posTPCNSigmaPr, float);
DECLARE_SOA_COLUMN(PosTPCNSigmaPi, posTPCNSigmaPi, float);
DECLARE_SOA_COLUMN(NegTPCNSigmaPr, negTPCNSigmaPr, float);
DECLARE_SOA_COLUMN(NegTPCNSigmaPi, negTPCNSigmaPi, float);

DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);

// Dynamic columns for jets (Px,Py,Pz):
DECLARE_SOA_DYNAMIC_COLUMN(JetPx, jetPx, //! Jet px
                           [](float jetPt, float jetPhi) -> float { return jetPt * std::cos(jetPhi); });
DECLARE_SOA_DYNAMIC_COLUMN(JetPy, jetPy, //! Jet py
                           [](float jetPt, float jetPhi) -> float { return jetPt * std::sin(jetPhi); });
DECLARE_SOA_DYNAMIC_COLUMN(JetPz, jetPz, //! Jet pz
                           [](float jetPt, float jetEta) -> float { return jetPt * std::sinh(jetEta); });
// Same for leading particles:
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePx, leadParticlePx, //! Leading particle px
                           [](float leadParticlePt, float leadParticlePhi) -> float { return leadParticlePt * std::cos(leadParticlePhi); });
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePy, leadParticlePy, //! Leading particle py
                           [](float leadParticlePt, float leadParticlePhi) -> float { return leadParticlePt * std::sin(leadParticlePhi); });
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePz, leadParticlePz, //! Leading particle pz
                           [](float leadParticlePt, float leadParticleEta) -> float { return leadParticlePt * std::sinh(leadParticleEta); });
} // namespace lambdajetpol

DECLARE_SOA_TABLE(RingCollisions, "AOD", "RINGCOLLISIONS",
                  o2::soa::Index<>,               // self-index: auto-assigned row number
                  lambdajetpol::CentFT0M,
                  lambdajetpol::CentFT0C,
                  lambdajetpol::CentFV0A);

namespace lambdajetpol
{
DECLARE_SOA_INDEX_COLUMN(RingCollision, ringCollision); // Declare index after table is available
} // namespace lambdajetpol

DECLARE_SOA_TABLE(RingJets, "AOD", "RINGJETS",
                  lambdajetpol::RingCollisionId,  // relational index -> RingCollisions
                  lambdajetpol::JetPt,
                  lambdajetpol::JetEta,
                  lambdajetpol::JetPhi,
                  lambdajetpol::JetNConstituents,
                  // Dynamic columns (explicitly bound to their static inputs):
                  lambdajetpol::JetPx<lambdajetpol::JetPt, lambdajetpol::JetPhi>,
                  lambdajetpol::JetPy<lambdajetpol::JetPt, lambdajetpol::JetPhi>,
                  lambdajetpol::JetPz<lambdajetpol::JetPt, lambdajetpol::JetEta>);

DECLARE_SOA_TABLE(RingLeadP, "AOD", "RINGLEADP",
                  lambdajetpol::RingCollisionId,
                  lambdajetpol::LeadParticlePt,
                  lambdajetpol::LeadParticleEta,
                  lambdajetpol::LeadParticlePhi,
                  // Dynamic columns:
                  lambdajetpol::LeadParticlePx<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticlePhi>,
                  lambdajetpol::LeadParticlePy<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticlePhi>,
                  lambdajetpol::LeadParticlePz<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticleEta>);

DECLARE_SOA_TABLE(RingLaV0s, "AOD", "RINGLAV0S",
                  lambdajetpol::RingCollisionId,
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
                  lambdajetpol::NegPhi,
                  lambdajetpol::PosTPCNSigmaPr,
                  lambdajetpol::PosTPCNSigmaPi,
                  lambdajetpol::NegTPCNSigmaPr,
                  lambdajetpol::NegTPCNSigmaPi,
                  lambdajetpol::V0CosPA,
                  lambdajetpol::V0Radius,
                  lambdajetpol::DcaV0Daughters,
                  lambdajetpol::DcaPosToPV,
                  lambdajetpol::DcaNegToPV);

using RingCollision = RingCollisions::iterator; // Useful shorthand
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LAMBDAJETPOLARIZATIONIONS_H_
