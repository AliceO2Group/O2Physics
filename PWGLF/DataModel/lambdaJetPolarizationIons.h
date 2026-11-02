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
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{

namespace lambdajetpol
{
// Collision information:
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);
DECLARE_SOA_COLUMN(Zvtx, zvtx, float);
DECLARE_SOA_COLUMN(InteractionRate, interactionRate, float);
// DECLARE_SOA_COLUMN(MagField, magField, float);

// Jet (and jet proxies) information:
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
// DECLARE_SOA_COLUMN(JetNConstituents, jetNConstituents, uint16_t); // Currently removed from datamodel.
// Other variables can better reveal jet quenching and help identify good selection criteria for quenched jets proxies

DECLARE_SOA_COLUMN(LeadParticlePt, leadParticlePt, float);
DECLARE_SOA_COLUMN(LeadParticleEta, leadParticleEta, float);
DECLARE_SOA_COLUMN(LeadParticlePhi, leadParticlePhi, float);

// V0 information:
DECLARE_SOA_COLUMN(V0Pt, v0Pt, float);
DECLARE_SOA_COLUMN(V0Eta, v0Eta, float);
DECLARE_SOA_COLUMN(V0Phi, v0Phi, float);

DECLARE_SOA_COLUMN(IsLambda, isLambda, bool); //! 0: antiLambda, 1: Lambda. There are no ambiguous candidates stored (those that pass both Lambda-specific and AntiLambda-specific checks)
// DECLARE_SOA_COLUMN(IsAntiLambda, isAntiLambda, bool);
DECLARE_SOA_COLUMN(MassV0, massV0, float);
// DECLARE_SOA_COLUMN(MassLambda, massLambda, float);
// DECLARE_SOA_COLUMN(MassAntiLambda, massAntiLambda, float);

DECLARE_SOA_COLUMN(PosPt, posPt, float); // Could consider rewriting this as proton/pion-like Pt, as in TPCNSigma, instead of Pos/Neg
DECLARE_SOA_COLUMN(PosEta, posEta, float);
DECLARE_SOA_COLUMN(PosPhi, posPhi, float);
DECLARE_SOA_COLUMN(NegPt, negPt, float);
DECLARE_SOA_COLUMN(NegEta, negEta, float);
DECLARE_SOA_COLUMN(NegPhi, negPhi, float);

// DECLARE_SOA_COLUMN(PosTPCNSigmaPr, posTPCNSigmaPr, float);
// DECLARE_SOA_COLUMN(PosTPCNSigmaPi, posTPCNSigmaPi, float);
// DECLARE_SOA_COLUMN(NegTPCNSigmaPr, negTPCNSigmaPr, float);
// DECLARE_SOA_COLUMN(NegTPCNSigmaPi, negTPCNSigmaPi, float);
// TPC Nsigma variables now stored only for non-ambiguous candidates (those that pass only the Lambda-specific or antiLambda-specific tests):
// (PrLike: proton-like daughter of the V0, either a proton (Lambda) or an antiproton (AntiLambda). PiLike is the pion daughter)
DECLARE_SOA_COLUMN(RoundPrLikeTPCNSigma, roundPrLikeTPCNSigma, int16_t); //! Stores nsigma with a precision of 0.01 nsigma. Following LFStrangenessPIDTables.h
DECLARE_SOA_COLUMN(RoundPiLikeTPCNSigma, roundPiLikeTPCNSigma, int16_t);

DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);

// Dynamic columns for jets (Px,Py,Pz):
DECLARE_SOA_DYNAMIC_COLUMN(JetPx, jetPx,
                           [](float jetPt, float jetPhi) -> float { return jetPt * std::cos(jetPhi); });
DECLARE_SOA_DYNAMIC_COLUMN(JetPy, jetPy,
                           [](float jetPt, float jetPhi) -> float { return jetPt * std::sin(jetPhi); });
DECLARE_SOA_DYNAMIC_COLUMN(JetPz, jetPz,
                           [](float jetPt, float jetEta) -> float { return jetPt * std::sinh(jetEta); });
// Same for leading particles:
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePx, leadParticlePx,
                           [](float leadParticlePt, float leadParticlePhi) -> float { return leadParticlePt * std::cos(leadParticlePhi); });
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePy, leadParticlePy,
                           [](float leadParticlePt, float leadParticlePhi) -> float { return leadParticlePt * std::sin(leadParticlePhi); });
DECLARE_SOA_DYNAMIC_COLUMN(LeadParticlePz, leadParticlePz,
                           [](float leadParticlePt, float leadParticleEta) -> float { return leadParticlePt * std::sinh(leadParticleEta); });

// Dynamic columns for retrieving rounded Nsigma columns:
DECLARE_SOA_DYNAMIC_COLUMN(PrLikeTPCNSigma, prLikeTPCNSigma,
                           [](int16_t packedValue) -> float { return packedValue / 100.f; });
DECLARE_SOA_DYNAMIC_COLUMN(PiLikeTPCNSigma, piLikeTPCNSigma,
                           [](int16_t packedValue) -> float { return packedValue / 100.f; });
} // namespace lambdajetpol

DECLARE_SOA_TABLE(RingCollisions, "AOD", "RINGCOLLISION",
                  o2::soa::Index<>, //! self-index: auto-assigned row number
                  lambdajetpol::CentFT0M,
                  lambdajetpol::CentFT0C,
                  lambdajetpol::CentFV0A,
                  lambdajetpol::Zvtx,
                  lambdajetpol::InteractionRate);

namespace lambdajetpol
{
DECLARE_SOA_INDEX_COLUMN(RingCollision, ringCollision); //! Declare index after table is available
} // namespace lambdajetpol

DECLARE_SOA_TABLE(RingJets, "AOD", "RINGJET",
                  lambdajetpol::RingCollisionId, //! relational index -> RingCollisions
                  lambdajetpol::JetPt,
                  lambdajetpol::JetEta,
                  lambdajetpol::JetPhi,
                  // lambdajetpol::JetNConstituents,
                  // Dynamic columns (explicitly bound to their static inputs):
                  lambdajetpol::JetPx<lambdajetpol::JetPt, lambdajetpol::JetPhi>,
                  lambdajetpol::JetPy<lambdajetpol::JetPt, lambdajetpol::JetPhi>,
                  lambdajetpol::JetPz<lambdajetpol::JetPt, lambdajetpol::JetEta>);

DECLARE_SOA_TABLE(RingLeadPs, "AOD", "RINGLEADP",
                  lambdajetpol::RingCollisionId,
                  lambdajetpol::LeadParticlePt,
                  lambdajetpol::LeadParticleEta,
                  lambdajetpol::LeadParticlePhi,
                  // Dynamic columns:
                  lambdajetpol::LeadParticlePx<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticlePhi>,
                  lambdajetpol::LeadParticlePy<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticlePhi>,
                  lambdajetpol::LeadParticlePz<lambdajetpol::LeadParticlePt, lambdajetpol::LeadParticleEta>);

DECLARE_SOA_TABLE(RingLaV0s, "AOD", "RINGLAV0",
                  lambdajetpol::RingCollisionId,
                  lambdajetpol::V0Pt,
                  lambdajetpol::V0Eta,
                  lambdajetpol::V0Phi,
                  lambdajetpol::IsLambda,
                  lambdajetpol::MassV0,
                  lambdajetpol::PosPt,
                  lambdajetpol::PosEta,
                  lambdajetpol::PosPhi,
                  lambdajetpol::NegPt,
                  lambdajetpol::NegEta,
                  lambdajetpol::NegPhi,
                  lambdajetpol::RoundPrLikeTPCNSigma,
                  lambdajetpol::RoundPiLikeTPCNSigma,
                  lambdajetpol::V0CosPA,
                  lambdajetpol::V0Radius,
                  lambdajetpol::DcaV0Daughters,
                  lambdajetpol::DcaPosToPV,
                  lambdajetpol::DcaNegToPV,
                  // Dynamic columns:
                  lambdajetpol::PrLikeTPCNSigma<lambdajetpol::RoundPrLikeTPCNSigma>,
                  lambdajetpol::PiLikeTPCNSigma<lambdajetpol::RoundPiLikeTPCNSigma>);

using RingCollision = RingCollisions::iterator; // Useful shorthand
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LAMBDAJETPOLARIZATIONIONS_H_
