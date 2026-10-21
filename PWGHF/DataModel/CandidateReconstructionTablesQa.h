// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CandidateReconstructionTablesQa.h
/// \brief QA candidate tables for Xic0/Omegac0 cascade decays
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston

#ifndef PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_
#define PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cstdint>

namespace o2::aod
{
namespace hf_cand_casc_bach
{
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
// indices to cascade and bachelor track
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);
DECLARE_SOA_INDEX_COLUMN_FULL(CharmBachelor, charmBachelor, int, Tracks, "_charmBachelor");
// sign of the decay, FIXME: probably can be removed, as it is redundant with the sign of the bachelor track
// DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t);
// vertices of decay chain
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float);
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float);
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float);
DECLARE_SOA_COLUMN(XCascadeVertex, xCascadeVertex, float);
DECLARE_SOA_COLUMN(YCascadeVertex, yCascadeVertex, float);
DECLARE_SOA_COLUMN(ZCascadeVertex, zCascadeVertex, float);
DECLARE_SOA_COLUMN(XV0Vertex, xV0Vertex, float);
DECLARE_SOA_COLUMN(YV0Vertex, yV0Vertex, float);
DECLARE_SOA_COLUMN(ZV0Vertex, zV0Vertex, float);
// masses, FIXME: check if they are necessary
// DECLARE_SOA_COLUMN(Mass, mass, float);
// DECLARE_SOA_COLUMN(MassCascade, massCascade, float);
// DECLARE_SOA_COLUMN(MassV0, massV0, float);
// momenta
DECLARE_SOA_COLUMN(PxCascade, pxCascade, float);
DECLARE_SOA_COLUMN(PyCascade, pyCascade, float);
DECLARE_SOA_COLUMN(PzCascade, pzCascade, float);
DECLARE_SOA_COLUMN(PxCharmBachelor, pxCharmBachelor, float);
DECLARE_SOA_COLUMN(PyCharmBachelor, pyCharmBachelor, float);
DECLARE_SOA_COLUMN(PzCharmBachelor, pzCharmBachelor, float);
DECLARE_SOA_COLUMN(PxV0, pxV0, float);
DECLARE_SOA_COLUMN(PyV0, pyV0, float);
DECLARE_SOA_COLUMN(PzV0, pzV0, float);
DECLARE_SOA_COLUMN(PxCascadeBachelor, pxCascadeBachelor, float);
DECLARE_SOA_COLUMN(PyCascadeBachelor, pyCascadeBachelor, float);
DECLARE_SOA_COLUMN(PzCascadeBachelor, pzCascadeBachelor, float);
DECLARE_SOA_COLUMN(PxV0PosDaughter, pxV0PosDaughter, float);
DECLARE_SOA_COLUMN(PyV0PosDaughter, pyV0PosDaughter, float);
DECLARE_SOA_COLUMN(PzV0PosDaughter, pzV0PosDaughter, float);
DECLARE_SOA_COLUMN(PxV0NegDaughter, pxV0NegDaughter, float);
DECLARE_SOA_COLUMN(PyV0NegDaughter, pyV0NegDaughter, float);
DECLARE_SOA_COLUMN(PzV0NegDaughter, pzV0NegDaughter, float);
// information from fitter
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float);
DECLARE_SOA_COLUMN(ImpactParameterCascadeXY, impactParameterCascadeXY, float);
DECLARE_SOA_COLUMN(ImpactParameterCascadeZ, impactParameterCascadeZ, float);
DECLARE_SOA_COLUMN(ErrorImpactParameterCascadeXY, errorImpactParameterCascadeXY, float);
DECLARE_SOA_COLUMN(ErrorImpactParameterCascadeZ, errorImpactParameterCascadeZ, float);
DECLARE_SOA_COLUMN(ImpactParameterCharmBachelorXY, impactParameterCharmBachelorXY, float);
DECLARE_SOA_COLUMN(ImpactParameterCharmBachelorZ, impactParameterCharmBachelorZ, float);
DECLARE_SOA_COLUMN(ErrorImpactParameterCharmBachelorXY, errorImpactParameterCharmBachelorXY, float);
DECLARE_SOA_COLUMN(ErrorImpactParameterCharmBachelorZ, errorImpactParameterCharmBachelorZ, float);

DECLARE_SOA_COLUMN(DcaXYCascadeBachelorToPv, dcaXYCascadeBachelorToPv, float);
DECLARE_SOA_COLUMN(DcaZCascadeBachelorToPv, dcaZCascadeBachelorToPv, float);
DECLARE_SOA_COLUMN(DcaXYV0PosDaughterToPv, dcaXYV0PosDaughterToPv, float);
DECLARE_SOA_COLUMN(DcaZV0PosDaughterToPv, dcaZV0PosDaughterToPv, float);
DECLARE_SOA_COLUMN(DcaXYV0NegDaughterToPv, dcaXYV0NegDaughterToPv, float);
DECLARE_SOA_COLUMN(DcaZV0NegDaughterToPv, dcaZV0NegDaughterToPv, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaCascadeDaughters, dcaCascadeDaughters, float);
DECLARE_SOA_COLUMN(DcaDaughters, dcaDaughters, float);

DECLARE_SOA_EXPRESSION_COLUMN(Px, px, float,
                              1.f * aod::hf_cand_casc_bach::pxCascade + 1.f * aod::hf_cand_casc_bach::pxCharmBachelor);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, float,
                              1.f * aod::hf_cand_casc_bach::pyCascade + 1.f * aod::hf_cand_casc_bach::pyCharmBachelor);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, float,
                              1.f * aod::hf_cand_casc_bach::pzCascade + 1.f * aod::hf_cand_casc_bach::pzCharmBachelor);

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Transverse momentum of the candidate
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Momentum of the candidate
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Pseudorapidity of the candidate
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Azimuthal angle of the candidate
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Rapidity for the supplied mass
                           [](float px, float py, float pz, double mass) -> float { return RecoDecay::y(std::array{px, py, pz}, mass); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //! Energy for the supplied mass
                           [](float px, float py, float pz, double mass) -> float { return RecoDecay::e(px, py, pz, mass); });
DECLARE_SOA_DYNAMIC_COLUMN(PtCascade, ptCascade, //! Transverse momentum of the cascade
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PCascade, pCascade, //! Momentum of the cascade
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaCascade, etaCascade, //! Pseudorapidity of the cascade
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiCascade, phiCascade, //! Azimuthal angle of the cascade
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0, ptV0, //! Transverse momentum of the V0
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PV0, pV0, //! Momentum of the V0
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaV0, etaV0, //! Pseudorapidity of the V0
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiV0, phiV0, //! Azimuthal angle of the V0
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaCharmBachelor, etaCharmBachelor, //! Pseudorapidity of the charm bachelor track
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtCharmBachelor, ptCharmBachelor, //! Transverse momentum of the charm bachelor track
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PCharmBachelor, pCharmBachelor, //! Momentum of the charm bachelor track
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiCharmBachelor, phiCharmBachelor, //! Azimuthal angle of the charm bachelor track
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaCascadeBachelor, etaCascadeBachelor, //! Pseudorapidity of the cascade bachelor track
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtCascadeBachelor, ptCascadeBachelor, //! Transverse momentum of the cascade bachelor track
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PCascadeBachelor, pCascadeBachelor, //! Momentum of the cascade bachelor track
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiCascadeBachelor, phiCascadeBachelor, //! Azimuthal angle of the cascade bachelor track
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaV0PosDaughter, etaV0PosDaughter, //! Pseudorapidity of the positive V0 daughter
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0PosDaughter, ptV0PosDaughter, //! Transverse momentum of the positive V0 daughter
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PV0PosDaughter, pV0PosDaughter, //! Momentum of the positive V0 daughter
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiV0PosDaughter, phiV0PosDaughter, //! Azimuthal angle of the positive V0 daughter
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(EtaV0NegDaughter, etaV0NegDaughter, //! Pseudorapidity of the negative V0 daughter
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(PtV0NegDaughter, ptV0NegDaughter, //! Transverse momentum of the negative V0 daughter
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PV0NegDaughter, pV0NegDaughter, //! Momentum of the negative V0 daughter
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PhiV0NegDaughter, phiV0NegDaughter, //! Azimuthal angle of the negative V0 daughter
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });

DECLARE_SOA_DYNAMIC_COLUMN(MV0, mV0, //! Invariant mass for supplied daughter-mass hypotheses
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& masses) -> float {
                             return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, masses);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MCascade, mCascade, //! Invariant mass for supplied daughter-mass hypotheses
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& masses) -> float {
                             return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, masses);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! Invariant mass for supplied daughter-mass hypotheses
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& masses) -> float {
                             return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, masses);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosThetaStarCharmBachelor, cosThetaStarCharmBachelor, //! Helicity angle for supplied daughter and parent mass hypotheses
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& masses, double parentMass, int prong) -> float {
                             return RecoDecay::cosThetaStar(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, masses, parentMass, prong);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //! Candidate decay length from the PV
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, //! Candidate transverse decay length from the PV
                           [](float pvx, float pvy, float svx, float svy) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, //! Candidate decay length divided by its uncertainty
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float error) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}) / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, //! Candidate transverse decay length divided by its uncertainty
                           [](float pvx, float pvy, float svx, float svy, float error) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy}) / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthCascade, decayLengthCascade, //! Cascade decay length from the parent vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYCascade, decayLengthXYCascade, //! Cascade transverse decay length from the parent vertex
                           [](float pvx, float pvy, float svx, float svy) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthV0, decayLengthV0, //! V0 decay length from the cascade vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYV0, decayLengthXYV0, //! V0 transverse decay length from the cascade vertex
                           [](float pvx, float pvy, float svx, float svy) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterCascadeXYNormalised, impactParameterCascadeXYNormalised, //! Cascade transverse impact parameter divided by its uncertainty
                           [](float impactParameter, float error) -> float {
                             return impactParameter / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterCascadeZNormalised, impactParameterCascadeZNormalised, //! Cascade longitudinal impact parameter divided by its uncertainty
                           [](float impactParameter, float error) -> float {
                             return impactParameter / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterCharmBachelorXYNormalised, impactParameterCharmBachelorXYNormalised, //! Charm-bachelor transverse impact parameter divided by its uncertainty
                           [](float impactParameter, float error) -> float {
                             return impactParameter / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterCharmBachelorZNormalised, impactParameterCharmBachelorZNormalised, //! Charm-bachelor longitudinal impact parameter divided by its uncertainty
                           [](float impactParameter, float error) -> float {
                             return impactParameter / error;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPA, cosPA, //! Candidate pointing angle from the PV
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz) -> float {
                             return RecoDecay::cpa(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}, std::array{mx, my, mz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAXY, cosPAXY, //! Candidate transverse pointing angle from the PV
                           [](float pvx, float pvy, float svx, float svy, float mx, float my) -> float {
                             return RecoDecay::cpaXY(std::array{pvx, pvy}, std::array{svx, svy}, std::array{mx, my});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPACascadeToPv, cosPACascadeToPv, //! Cascade pointing angle from the PV
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz) -> float {
                             return RecoDecay::cpa(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}, std::array{mx, my, mz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAXYCascadeToPv, cosPAXYCascadeToPv, //! Cascade transverse pointing angle from the PV
                           [](float pvx, float pvy, float svx, float svy, float mx, float my) -> float {
                             return RecoDecay::cpaXY(std::array{pvx, pvy}, std::array{svx, svy}, std::array{mx, my});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAV0ToPv, cosPAV0ToPv, //! V0 pointing angle from the PV
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz) -> float {
                             return RecoDecay::cpa(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}, std::array{mx, my, mz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAXYV0ToPv, cosPAXYV0ToPv, //! V0 transverse pointing angle from the PV
                           [](float pvx, float pvy, float svx, float svy, float mx, float my) -> float {
                             return RecoDecay::cpaXY(std::array{pvx, pvy}, std::array{svx, svy}, std::array{mx, my});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAV0ToCascade, cosPAV0ToCascade, //! V0 pointing angle from the cascade vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz) -> float {
                             return RecoDecay::cpa(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}, std::array{mx, my, mz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAXYV0ToCascade, cosPAXYV0ToCascade, //! V0 transverse pointing angle from the cascade vertex
                           [](float pvx, float pvy, float svx, float svy, float mx, float my) -> float {
                             return RecoDecay::cpaXY(std::array{pvx, pvy}, std::array{svx, svy}, std::array{mx, my});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPACascadeToParent, cosPACascadeToParent, //! Cascade pointing angle from the parent vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz) -> float {
                             return RecoDecay::cpa(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}, std::array{mx, my, mz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CosPAXYCascadeToParent, cosPAXYCascadeToParent, //! Cascade transverse pointing angle from the parent vertex
                           [](float pvx, float pvy, float svx, float svy, float mx, float my) -> float {
                             return RecoDecay::cpaXY(std::array{pvx, pvy}, std::array{svx, svy}, std::array{mx, my});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Ct, ct, //! Parent proper decay length for the supplied parent mass
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz, double mass) -> float {
                             return RecoDecay::ct(std::array{mx, my, mz}, RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}), mass);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CtCascade, ctCascade, //! Cascade proper decay length for the supplied cascade mass
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz, double mass) -> float {
                             return RecoDecay::ct(std::array{mx, my, mz}, RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}), mass);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(CtV0, ctV0, //! V0 proper decay length for the supplied V0 mass
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz, float mx, float my, float mz, double mass) -> float {
                             return RecoDecay::ct(std::array{mx, my, mz}, RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz}), mass);
                           });

// KFParticle specific columns
DECLARE_SOA_COLUMN(MassErrorV0, massErrorV0, float);
DECLARE_SOA_COLUMN(MassCascadeAlternative, massCascadeAlternative, float);
DECLARE_SOA_COLUMN(MassErrorCascade, massErrorCascade, float);
DECLARE_SOA_COLUMN(MassErrorCascadeAlternative, massErrorCascadeAlternative, float);
DECLARE_SOA_COLUMN(MassErrorParent, massErrorParent, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCascade, chi2GeoCascade, float);
DECLARE_SOA_COLUMN(Chi2GeoParent, chi2GeoParent, float);
DECLARE_SOA_COLUMN(NdfGeoV0, ndfGeoV0, int32_t);
DECLARE_SOA_COLUMN(NdfGeoCascade, ndfGeoCascade, int32_t);
DECLARE_SOA_COLUMN(NdfGeoParent, ndfGeoParent, int32_t);
DECLARE_SOA_COLUMN(Chi2MassV0, chi2MassV0, float);
DECLARE_SOA_COLUMN(Chi2MassCascade, chi2MassCascade, float);
DECLARE_SOA_COLUMN(NdfMassV0, ndfMassV0, int32_t);
DECLARE_SOA_COLUMN(NdfMassCascade, ndfMassCascade, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(NdfTopoV0ToPv, ndfTopoV0ToPv, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCascadeToPv, chi2TopoCascadeToPv, float);
DECLARE_SOA_COLUMN(NdfTopoCascadeToPv, ndfTopoCascadeToPv, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCharmBachelorToPv, chi2TopoCharmBachelorToPv, float);
DECLARE_SOA_COLUMN(NdfTopoCharmBachelorToPv, ndfTopoCharmBachelorToPv, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoParentToPv, chi2TopoParentToPv, float);
DECLARE_SOA_COLUMN(NdfTopoParentToPv, ndfTopoParentToPv, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCascade, chi2TopoV0ToCascade, float);
DECLARE_SOA_COLUMN(NdfTopoV0ToCascade, ndfTopoV0ToCascade, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCascadeToParent, chi2TopoCascadeToParent, float);
DECLARE_SOA_COLUMN(NdfTopoCascadeToParent, ndfTopoCascadeToParent, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCascadeBachelorToCascade, chi2TopoCascadeBachelorToCascade, float);
DECLARE_SOA_COLUMN(NdfTopoCascadeBachelorToCascade, ndfTopoCascadeBachelorToCascade, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCharmBachelorToParent, chi2TopoCharmBachelorToParent, float);
DECLARE_SOA_COLUMN(NdfTopoCharmBachelorToParent, ndfTopoCharmBachelorToParent, int32_t);
DECLARE_SOA_COLUMN(DecayLengthSignificanceV0, decayLengthSignificanceV0, float);
DECLARE_SOA_COLUMN(DecayLengthSignificanceCascade, decayLengthSignificanceCascade, float);
DECLARE_SOA_COLUMN(DecayLengthSignificanceParent, decayLengthSignificanceParent, float);
DECLARE_SOA_COLUMN(DcaXYCascadeToPv, dcaXYCascadeToPv, float);
DECLARE_SOA_COLUMN(DcaXYCharmBachelorToPv, dcaXYCharmBachelorToPv, float);
DECLARE_SOA_COLUMN(DeviationCharmBachelorToPv, deviationCharmBachelorToPv, float);

DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdfV0, chi2GeoOverNdfV0, //! Geometric-fit chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdfCascade, chi2GeoOverNdfCascade, //! Geometric-fit chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdfParent, chi2GeoOverNdfParent, //! Geometric-fit chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2MassOverNdfV0, chi2MassOverNdfV0, //! Mass-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2MassOverNdfCascade, chi2MassOverNdfCascade, //! Mass-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfV0ToPv, chi2TopoOverNdfV0ToPv, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCascadeToPv, chi2TopoOverNdfCascadeToPv, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCharmBachelorToPv, chi2TopoOverNdfCharmBachelorToPv, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfParentToPv, chi2TopoOverNdfParentToPv, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfV0ToCascade, chi2TopoOverNdfV0ToCascade, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCascadeToParent, chi2TopoOverNdfCascadeToParent, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCascadeBachelorToCascade, chi2TopoOverNdfCascadeBachelorToCascade, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCharmBachelorToParent, chi2TopoOverNdfCharmBachelorToParent, //! Topological-constraint chi2 per NDF
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });

} // namespace hf_cand_casc_bach

DECLARE_SOA_TABLE(HfCandCascBach, "AOD", "HFCANDCASCBACH", //! Common QA candidate table for XiPi, OmegaPi and OmegaK
                  o2::soa::Index<>,
                  hf_cand_casc_bach::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  // hf_cand_casc_bach::SignDecay,
                  hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex,
                  hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex,
                  hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::ZV0Vertex,
                  hf_cand_casc_bach::ErrorDecayLength, hf_cand_casc_bach::ErrorDecayLengthXY,
                  hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade,
                  hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor, hf_cand_casc_bach::PzCharmBachelor,
                  hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0,
                  hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor, hf_cand_casc_bach::PzCascadeBachelor,
                  hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter, hf_cand_casc_bach::PzV0PosDaughter,
                  hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter, hf_cand_casc_bach::PzV0NegDaughter,
                  // hf_cand_casc_bach::MassV0, hf_cand_casc_bach::MassCascade, hf_cand_casc_bach::Mass,
                  hf_cand_casc_bach::ImpactParameterCascadeXY, hf_cand_casc_bach::ImpactParameterCascadeZ,
                  hf_cand_casc_bach::ErrorImpactParameterCascadeXY, hf_cand_casc_bach::ErrorImpactParameterCascadeZ,
                  hf_cand_casc_bach::ImpactParameterCharmBachelorXY, hf_cand_casc_bach::ImpactParameterCharmBachelorZ,
                  hf_cand_casc_bach::ErrorImpactParameterCharmBachelorXY, hf_cand_casc_bach::ErrorImpactParameterCharmBachelorZ,
                  hf_cand_casc_bach::DcaXYV0PosDaughterToPv, hf_cand_casc_bach::DcaXYV0NegDaughterToPv, hf_cand_casc_bach::DcaXYCascadeBachelorToPv,
                  hf_cand_casc_bach::DcaZV0PosDaughterToPv, hf_cand_casc_bach::DcaZV0NegDaughterToPv, hf_cand_casc_bach::DcaZCascadeBachelorToPv,
                  hf_cand_casc_bach::DcaV0Daughters, hf_cand_casc_bach::DcaCascadeDaughters, hf_cand_casc_bach::DcaDaughters,
                  hf_cand_casc_bach::CascadeId, hf_cand_casc_bach::CharmBachelorId,
                  hf_track_index::HFflag,
                  /* dynamic and expression columns */
                  hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz,
                  hf_cand_casc_bach::Pt<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py>,
                  hf_cand_casc_bach::P<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::Eta<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::Phi<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py>,
                  hf_cand_casc_bach::Y<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::E<hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::PtCascade<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade>,
                  hf_cand_casc_bach::PCascade<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade>,
                  hf_cand_casc_bach::EtaCascade<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade>,
                  hf_cand_casc_bach::PhiCascade<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade>,
                  hf_cand_casc_bach::PtV0<hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0>,
                  hf_cand_casc_bach::PV0<hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0>,
                  hf_cand_casc_bach::EtaV0<hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0>,
                  hf_cand_casc_bach::PhiV0<hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0>,
                  hf_cand_casc_bach::EtaCharmBachelor<hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor, hf_cand_casc_bach::PzCharmBachelor>,
                  hf_cand_casc_bach::PtCharmBachelor<hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor>,
                  hf_cand_casc_bach::PCharmBachelor<hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor, hf_cand_casc_bach::PzCharmBachelor>,
                  hf_cand_casc_bach::PhiCharmBachelor<hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor>,
                  hf_cand_casc_bach::EtaCascadeBachelor<hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor, hf_cand_casc_bach::PzCascadeBachelor>,
                  hf_cand_casc_bach::PtCascadeBachelor<hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor>,
                  hf_cand_casc_bach::PCascadeBachelor<hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor, hf_cand_casc_bach::PzCascadeBachelor>,
                  hf_cand_casc_bach::PhiCascadeBachelor<hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor>,
                  hf_cand_casc_bach::EtaV0PosDaughter<hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter, hf_cand_casc_bach::PzV0PosDaughter>,
                  hf_cand_casc_bach::PtV0PosDaughter<hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter>,
                  hf_cand_casc_bach::PV0PosDaughter<hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter, hf_cand_casc_bach::PzV0PosDaughter>,
                  hf_cand_casc_bach::PhiV0PosDaughter<hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter>,
                  hf_cand_casc_bach::EtaV0NegDaughter<hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter, hf_cand_casc_bach::PzV0NegDaughter>,
                  hf_cand_casc_bach::PtV0NegDaughter<hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter>,
                  hf_cand_casc_bach::PV0NegDaughter<hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter, hf_cand_casc_bach::PzV0NegDaughter>,
                  hf_cand_casc_bach::PhiV0NegDaughter<hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter>,
                  hf_cand_casc_bach::MV0<hf_cand_casc_bach::PxV0PosDaughter, hf_cand_casc_bach::PyV0PosDaughter, hf_cand_casc_bach::PzV0PosDaughter, hf_cand_casc_bach::PxV0NegDaughter, hf_cand_casc_bach::PyV0NegDaughter, hf_cand_casc_bach::PzV0NegDaughter>,
                  hf_cand_casc_bach::MCascade<hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0, hf_cand_casc_bach::PxCascadeBachelor, hf_cand_casc_bach::PyCascadeBachelor, hf_cand_casc_bach::PzCascadeBachelor>,
                  hf_cand_casc_bach::M<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade, hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor, hf_cand_casc_bach::PzCharmBachelor>,
                  hf_cand_casc_bach::CosThetaStarCharmBachelor<hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade, hf_cand_casc_bach::PxCharmBachelor, hf_cand_casc_bach::PyCharmBachelor, hf_cand_casc_bach::PzCharmBachelor>,
                  hf_cand_casc_bach::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex>,
                  hf_cand_casc_bach::DecayLengthXY<collision::PosX, collision::PosY, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex>,
                  hf_cand_casc_bach::DecayLengthCascade<hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex>,
                  hf_cand_casc_bach::DecayLengthXYCascade<hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex>,
                  hf_cand_casc_bach::DecayLengthV0<hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::ZV0Vertex>,
                  hf_cand_casc_bach::DecayLengthXYV0<hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex>,
                  hf_cand_casc_bach::CosPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::CosPAXY<collision::PosX, collision::PosY, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::Px, hf_cand_casc_bach::Py>,
                  hf_cand_casc_bach::CosPACascadeToPv<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade>,
                  hf_cand_casc_bach::CosPAXYCascadeToPv<collision::PosX, collision::PosY, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade>,
                  hf_cand_casc_bach::CosPAV0ToPv<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::ZV0Vertex, hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0>,
                  hf_cand_casc_bach::CosPAXYV0ToPv<collision::PosX, collision::PosY, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0>,
                  hf_cand_casc_bach::CosPAV0ToCascade<hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::ZV0Vertex, hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0>,
                  hf_cand_casc_bach::CosPAXYV0ToCascade<hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0>,
                  hf_cand_casc_bach::CosPACascadeToParent<hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade>,
                  hf_cand_casc_bach::CosPAXYCascadeToParent<hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade>,
                  hf_cand_casc_bach::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::Px, hf_cand_casc_bach::Py, hf_cand_casc_bach::Pz>,
                  hf_cand_casc_bach::CtCascade<hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::PxCascade, hf_cand_casc_bach::PyCascade, hf_cand_casc_bach::PzCascade>,
                  hf_cand_casc_bach::CtV0<hf_cand_casc_bach::XCascadeVertex, hf_cand_casc_bach::YCascadeVertex, hf_cand_casc_bach::ZCascadeVertex, hf_cand_casc_bach::XV0Vertex, hf_cand_casc_bach::YV0Vertex, hf_cand_casc_bach::ZV0Vertex, hf_cand_casc_bach::PxV0, hf_cand_casc_bach::PyV0, hf_cand_casc_bach::PzV0>,
                  hf_cand_casc_bach::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ZSecondaryVertex, hf_cand_casc_bach::ErrorDecayLength>,
                  hf_cand_casc_bach::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand_casc_bach::XSecondaryVertex, hf_cand_casc_bach::YSecondaryVertex, hf_cand_casc_bach::ErrorDecayLengthXY>,
                  hf_cand_casc_bach::ImpactParameterCascadeXYNormalised<hf_cand_casc_bach::ImpactParameterCascadeXY, hf_cand_casc_bach::ErrorImpactParameterCascadeXY>,
                  hf_cand_casc_bach::ImpactParameterCascadeZNormalised<hf_cand_casc_bach::ImpactParameterCascadeZ, hf_cand_casc_bach::ErrorImpactParameterCascadeZ>,
                  hf_cand_casc_bach::ImpactParameterCharmBachelorXYNormalised<hf_cand_casc_bach::ImpactParameterCharmBachelorXY, hf_cand_casc_bach::ErrorImpactParameterCharmBachelorXY>,
                  hf_cand_casc_bach::ImpactParameterCharmBachelorZNormalised<hf_cand_casc_bach::ImpactParameterCharmBachelorZ, hf_cand_casc_bach::ErrorImpactParameterCharmBachelorZ>);

DECLARE_SOA_TABLE(HfCandCascBachKf, "AOD", "HFCANDCASCBACHKF", //! KF specific information for HfCandCascBach
                  hf_cand_casc_bach::MassErrorV0,
                  hf_cand_casc_bach::MassCascadeAlternative,
                  hf_cand_casc_bach::MassErrorCascade,
                  hf_cand_casc_bach::MassErrorCascadeAlternative, hf_cand_casc_bach::MassErrorParent,
                  hf_cand_casc_bach::Chi2GeoV0, hf_cand_casc_bach::NdfGeoV0,
                  hf_cand_casc_bach::Chi2GeoCascade, hf_cand_casc_bach::NdfGeoCascade,
                  hf_cand_casc_bach::Chi2GeoParent, hf_cand_casc_bach::NdfGeoParent,
                  hf_cand_casc_bach::Chi2MassV0, hf_cand_casc_bach::NdfMassV0,
                  hf_cand_casc_bach::Chi2MassCascade, hf_cand_casc_bach::NdfMassCascade,
                  hf_cand_casc_bach::Chi2TopoV0ToPv, hf_cand_casc_bach::NdfTopoV0ToPv,
                  hf_cand_casc_bach::Chi2TopoCascadeToPv, hf_cand_casc_bach::NdfTopoCascadeToPv,
                  hf_cand_casc_bach::Chi2TopoCharmBachelorToPv, hf_cand_casc_bach::NdfTopoCharmBachelorToPv,
                  hf_cand_casc_bach::Chi2TopoParentToPv, hf_cand_casc_bach::NdfTopoParentToPv,
                  hf_cand_casc_bach::Chi2TopoV0ToCascade, hf_cand_casc_bach::NdfTopoV0ToCascade,
                  hf_cand_casc_bach::Chi2TopoCascadeToParent, hf_cand_casc_bach::NdfTopoCascadeToParent,
                  hf_cand_casc_bach::Chi2TopoCascadeBachelorToCascade, hf_cand_casc_bach::NdfTopoCascadeBachelorToCascade,
                  hf_cand_casc_bach::Chi2TopoCharmBachelorToParent, hf_cand_casc_bach::NdfTopoCharmBachelorToParent,
                  hf_cand_casc_bach::DecayLengthSignificanceV0,
                  hf_cand_casc_bach::DecayLengthSignificanceCascade,
                  hf_cand_casc_bach::DecayLengthSignificanceParent,
                  hf_cand_casc_bach::DcaXYCascadeToPv, hf_cand_casc_bach::DcaXYCharmBachelorToPv,
                  hf_cand_casc_bach::DeviationCharmBachelorToPv,
                  /* derived ratios */
                  hf_cand_casc_bach::Chi2GeoOverNdfV0<hf_cand_casc_bach::Chi2GeoV0, hf_cand_casc_bach::NdfGeoV0>,
                  hf_cand_casc_bach::Chi2GeoOverNdfCascade<hf_cand_casc_bach::Chi2GeoCascade, hf_cand_casc_bach::NdfGeoCascade>,
                  hf_cand_casc_bach::Chi2GeoOverNdfParent<hf_cand_casc_bach::Chi2GeoParent, hf_cand_casc_bach::NdfGeoParent>,
                  hf_cand_casc_bach::Chi2MassOverNdfV0<hf_cand_casc_bach::Chi2MassV0, hf_cand_casc_bach::NdfMassV0>,
                  hf_cand_casc_bach::Chi2MassOverNdfCascade<hf_cand_casc_bach::Chi2MassCascade, hf_cand_casc_bach::NdfMassCascade>,
                  hf_cand_casc_bach::Chi2TopoOverNdfV0ToPv<hf_cand_casc_bach::Chi2TopoV0ToPv, hf_cand_casc_bach::NdfTopoV0ToPv>,
                  hf_cand_casc_bach::Chi2TopoOverNdfCascadeToPv<hf_cand_casc_bach::Chi2TopoCascadeToPv, hf_cand_casc_bach::NdfTopoCascadeToPv>,
                  hf_cand_casc_bach::Chi2TopoOverNdfCharmBachelorToPv<hf_cand_casc_bach::Chi2TopoCharmBachelorToPv, hf_cand_casc_bach::NdfTopoCharmBachelorToPv>,
                  hf_cand_casc_bach::Chi2TopoOverNdfParentToPv<hf_cand_casc_bach::Chi2TopoParentToPv, hf_cand_casc_bach::NdfTopoParentToPv>,
                  hf_cand_casc_bach::Chi2TopoOverNdfV0ToCascade<hf_cand_casc_bach::Chi2TopoV0ToCascade, hf_cand_casc_bach::NdfTopoV0ToCascade>,
                  hf_cand_casc_bach::Chi2TopoOverNdfCascadeToParent<hf_cand_casc_bach::Chi2TopoCascadeToParent, hf_cand_casc_bach::NdfTopoCascadeToParent>,
                  hf_cand_casc_bach::Chi2TopoOverNdfCascadeBachelorToCascade<hf_cand_casc_bach::Chi2TopoCascadeBachelorToCascade, hf_cand_casc_bach::NdfTopoCascadeBachelorToCascade>,
                  hf_cand_casc_bach::Chi2TopoOverNdfCharmBachelorToParent<hf_cand_casc_bach::Chi2TopoCharmBachelorToParent, hf_cand_casc_bach::NdfTopoCharmBachelorToParent>);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandCascBachMCRec, "AOD", "HFCANDCASCBACHMCREC", //!
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::OriginMcRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::PtBhadMotherPart,
                  hf_cand_mc_flag::PdgBhadMotherPart,
                  hf_cand_mc_flag::NTracksDecayed,
                  hf_cand_mc_flag::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandCascBachMCGen, "AOD", "HFCANDCASCBACHMCGEN", //!
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_cand_mc_flag::OriginMcGen,
                  hf_cand_mc_flag::FlagMcDecayChanGen,
                  hf_cand_mc_flag::IdxBhadMotherPart);

} // namespace o2::aod

struct HfHelperQa {
  template <typename T>
  static auto invMassLambda(const T& candidate)
  {
    return candidate.mV0(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassAntiLambda(const T& candidate)
  {
    return candidate.mV0(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassProton});
  }

  template <typename T>
  static auto invMassXi(const T& candidate)
  {
    return candidate.mCascade(std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmega(const T& candidate)
  {
    return candidate.mCascade(std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto invMassXiczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToOmegaK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto yXiczero(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto eXiczero(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto ctXiczero(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto cosThetaStarXiczeroToXiPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassXiC0, 1);
  }

  template <typename T>
  static auto yOmegaczero(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto eOmegaczero(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto ctOmegaczero(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToXiPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassOmegaC0, 1);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassOmegaC0, 1);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToOmegaK(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassKPlus}, o2::constants::physics::MassOmegaC0, 1);
  }

  template <typename T>
  static auto ctXi(const T& candidate)
  {
    return candidate.ctCascade(o2::constants::physics::MassXiMinus);
  }

  template <typename T>
  static auto ctOmega(const T& candidate)
  {
    return candidate.ctCascade(o2::constants::physics::MassOmegaMinus);
  }

  template <typename T>
  static auto ctLambda(const T& candidate)
  {
    return candidate.ctV0(o2::constants::physics::MassLambda0);
  }
};

namespace o2::hf_decay
{
using HfDecayChannel = int8_t;

/// @brief XiC0 and OmegaC0 candidates: main channels
namespace hf_cand_casc_bach
{
enum DecayChannelMain : HfDecayChannel {
  // XiC0
  Xic0ToXiPi = 1,   // Xi- π+ → Λ π- π+ → p π- π- π+
  Xic0ToOmegaK = 2, // Ω- K+ → Λ K- K+ → p π- K- K+
  // OmegaC0
  Omegac0ToOmegaPi = 3,         // Ω- π+ → Λ π- π+ → p π- π- π+
  Omegac0ToOmegaK = 4,          // Ω- K+ → Λ K- K+ → p π- K- K+
  Omegac0ToXiPi = 5,            // Xi- π+ → Λ π- π+ → p π- π- π+
  NChannelsMain = Omegac0ToXiPi // last channel
};
} // namespace hf_cand_casc_bach

} // namespace o2::hf_decay

#endif // PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_
