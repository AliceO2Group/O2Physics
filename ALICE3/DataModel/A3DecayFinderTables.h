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
/// \file   A3DecayFinderTables.h
/// \since  04/07/2023
/// \brief  Set of tables for ALICE 3 decay finder
///

#ifndef ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_
#define ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_

// O2 includes
#include "Common/Core/RecoDecay.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"

enum a3selectionBit : uint32_t { kDCAxy = 0,
                                 kInnerTOFPion,
                                 kInnerTOFKaon,
                                 kInnerTOFProton,
                                 kOuterTOFPion,
                                 kOuterTOFKaon,
                                 kOuterTOFProton,
                                 kRICHPion,
                                 kRICHKaon,
                                 kRICHProton,
                                 kTruePion,
                                 kTrueKaon,
                                 kTrueProton,
                                 kTruePiPlusFromD,
                                 kTrueKaPlusFromD,
                                 kTruePiMinusFromD,
                                 kTrueKaMinusFromD,
                                 kTruePiPlusFromLc,
                                 kTrueKaPlusFromLc,
                                 kTruePrPlusFromLc,
                                 kTruePiMinusFromLc,
                                 kTrueKaMinusFromLc,
                                 kTruePrMinusFromLc,
                                 kTrueXiFromXiC,
                                 kTruePiFromXiC,
                                 kTruePiFromXiCC };

namespace o2::aod
{

// general decay properties
namespace a3_hf_cand
{
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, double); //!
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, double); //!
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, double); //!
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex,  //!
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float); //! sum of (non-weighted) distances of the secondary vertex to its prongs
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //!
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong0, pt2Prong0, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng0, pVectorProng0, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameterY0, impactParameterY0, float);                   //!
DECLARE_SOA_COLUMN(ErrorImpactParameterY0, errorImpactParameterY0, float);         //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ0, impactParameterZ0, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ0, errorImpactParameterZ0, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised0, impactParameterZNormalised0, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong1, pt2Prong1, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng1, pVectorProng1, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameterY1, impactParameterY1, float);                   //!
DECLARE_SOA_COLUMN(ErrorImpactParameterY1, errorImpactParameterY1, float);         //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ1, impactParameterZ1, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ1, errorImpactParameterZ1, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised1, impactParameterZNormalised1, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float); //!
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float); //!
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng2, ptProng2, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2Prong2, pt2Prong2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PVectorProng2, pVectorProng2, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_COLUMN(ImpactParameterY2, impactParameterY2, float);                   //!
DECLARE_SOA_COLUMN(ErrorImpactParameterY2, errorImpactParameterY2, float);         //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });
DECLARE_SOA_COLUMN(ImpactParameterZ2, impactParameterZ2, float);                     //!
DECLARE_SOA_COLUMN(ErrorImpactParameterZ2, errorImpactParameterZ2, float);           //!
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterZNormalised2, impactParameterZNormalised2, //!
                           [](float dca, float err) -> float { return dca / err; });

/// prong PID nsigma
DECLARE_SOA_COLUMN(NSigTrkPi0, nSigTrkPi0, float);       //!
DECLARE_SOA_COLUMN(NSigTrkKa0, nSigTrkKa0, float);       //!
DECLARE_SOA_COLUMN(NSigTrkPr0, nSigTrkPr0, float);       //!
DECLARE_SOA_COLUMN(NSigTrkPi1, nSigTrkPi1, float);       //!
DECLARE_SOA_COLUMN(NSigTrkKa1, nSigTrkKa1, float);       //!
DECLARE_SOA_COLUMN(NSigTrkPr1, nSigTrkPr1, float);       //!
DECLARE_SOA_COLUMN(NSigTrkPi2, nSigTrkPi2, float);       //!
DECLARE_SOA_COLUMN(NSigTrkKa2, nSigTrkKa2, float);       //!
DECLARE_SOA_COLUMN(NSigTrkPr2, nSigTrkPr2, float);       //!
DECLARE_SOA_COLUMN(NSigRichPi0, nSigRichPi0, float);     //!
DECLARE_SOA_COLUMN(NSigRichKa0, nSigRichKa0, float);     //!
DECLARE_SOA_COLUMN(NSigRichPr0, nSigRichPr0, float);     //!
DECLARE_SOA_COLUMN(NSigRichPi1, nSigRichPi1, float);     //!
DECLARE_SOA_COLUMN(NSigRichKa1, nSigRichKa1, float);     //!
DECLARE_SOA_COLUMN(NSigRichPr1, nSigRichPr1, float);     //!
DECLARE_SOA_COLUMN(NSigRichPi2, nSigRichPi2, float);     //!
DECLARE_SOA_COLUMN(NSigRichKa2, nSigRichKa2, float);     //!
DECLARE_SOA_COLUMN(NSigRichPr2, nSigRichPr2, float);     //!
DECLARE_SOA_COLUMN(NSigInnTofPi0, nSigInnTofPi0, float); //!
DECLARE_SOA_COLUMN(NSigInnTofKa0, nSigInnTofKa0, float); //!
DECLARE_SOA_COLUMN(NSigInnTofPr0, nSigInnTofPr0, float); //!
DECLARE_SOA_COLUMN(NSigInnTofPi1, nSigInnTofPi1, float); //!
DECLARE_SOA_COLUMN(NSigInnTofKa1, nSigInnTofKa1, float); //!
DECLARE_SOA_COLUMN(NSigInnTofPr1, nSigInnTofPr1, float); //!
DECLARE_SOA_COLUMN(NSigInnTofPi2, nSigInnTofPi2, float); //!
DECLARE_SOA_COLUMN(NSigInnTofKa2, nSigInnTofKa2, float); //!
DECLARE_SOA_COLUMN(NSigInnTofPr2, nSigInnTofPr2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi0, nSigOutTofPi0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa0, nSigOutTofKa0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr0, nSigOutTofPr0, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi1, nSigOutTofPi1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa1, nSigOutTofKa1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr1, nSigOutTofPr1, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPi2, nSigOutTofPi2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofKa2, nSigOutTofKa2, float); //!
DECLARE_SOA_COLUMN(NSigOutTofPr2, nSigOutTofPr2, float); //!

// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt2, pt2, //!
                           [](float px, float py) -> float { return RecoDecay::pt2(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(P2, p2, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::p2(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PVector, pVector, //!
                           [](float px, float py, float pz) -> std::array<float, 3> { return std::array{px, py, pz}; });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //!
                           [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::y(std::array{px, py, pz}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(E2, e2, //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e2(px, py, pz, m); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXY, decayLengthXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthNormalised, decayLengthNormalised, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float err) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}) / err; });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float err) -> float { return RecoDecay::distanceXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}) / err; });
DECLARE_SOA_COLUMN(ErrorDecayLength, errorDecayLength, float);     //!
DECLARE_SOA_COLUMN(ErrorDecayLengthXY, errorDecayLengthXY, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(Cpa, cpa,                               //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(CpaXY, cpaXY, //!
                           [](float xVtxP, float yVtxP, float xVtxS, float yVtxS, float px, float py) -> float { return RecoDecay::cpaXY(std::array{xVtxP, yVtxP}, std::array{xVtxS, yVtxS}, std::array{px, py}); });
DECLARE_SOA_DYNAMIC_COLUMN(Ct, ct, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz, double m) -> float { return RecoDecay::ct(std::array{px, py, pz}, RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}), m); });
DECLARE_SOA_DYNAMIC_COLUMN(ImpactParameterXY, impactParameterXY, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::impParXY(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });

// ML scores columns
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float); //!
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float); //!
DECLARE_SOA_COLUMN(MlScore2, mlScore2, float); //!
} // namespace a3_hf_cand

#define HFCAND_COLUMNS                                                                                                                                                                                            \
  a3_hf_cand::CollisionId,                                                                                                                                                                                        \
    collision::PosX, collision::PosY, collision::PosZ,                                                                                                                                                            \
    a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ZSecondaryVertex,                                                                                                                     \
    a3_hf_cand::ErrorDecayLength, a3_hf_cand::ErrorDecayLengthXY,                                                                                                                                                 \
    a3_hf_cand::Chi2PCA,                                                                                                                                                                                          \
    /* dynamic columns */ a3_hf_cand::RSecondaryVertex<a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex>,                                                                                               \
    a3_hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ZSecondaryVertex>,                                         \
    a3_hf_cand::DecayLengthXY<collision::PosX, collision::PosY, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex>,                                                                                      \
    a3_hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ZSecondaryVertex, a3_hf_cand::ErrorDecayLength>, \
    a3_hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ErrorDecayLengthXY>

// general columns
#define HFPRONG0_COLUMNS                                                                                     \
  a3_hf_cand::ImpactParameterNormalised0<a3_hf_cand::ImpactParameterY0, a3_hf_cand::ErrorImpactParameterY0>, \
    a3_hf_cand::PtProng0<a3_hf_cand::PxProng0, a3_hf_cand::PyProng0>,                                        \
    a3_hf_cand::Pt2Prong0<a3_hf_cand::PxProng0, a3_hf_cand::PyProng0>,                                       \
    a3_hf_cand::PVectorProng0<a3_hf_cand::PxProng0, a3_hf_cand::PyProng0, a3_hf_cand::PzProng0>

#define HFPRONG1_COLUMNS                                                                                     \
  a3_hf_cand::ImpactParameterNormalised1<a3_hf_cand::ImpactParameterY1, a3_hf_cand::ErrorImpactParameterY1>, \
    a3_hf_cand::PtProng1<a3_hf_cand::PxProng1, a3_hf_cand::PyProng1>,                                        \
    a3_hf_cand::Pt2Prong1<a3_hf_cand::PxProng1, a3_hf_cand::PyProng1>,                                       \
    a3_hf_cand::PVectorProng1<a3_hf_cand::PxProng1, a3_hf_cand::PyProng1, a3_hf_cand::PzProng1>

#define HFPRONG2_COLUMNS                                                                                     \
  a3_hf_cand::ImpactParameterNormalised2<a3_hf_cand::ImpactParameterY2, a3_hf_cand::ErrorImpactParameterY2>, \
    a3_hf_cand::PtProng2<a3_hf_cand::PxProng2, a3_hf_cand::PyProng2>,                                        \
    a3_hf_cand::Pt2Prong2<a3_hf_cand::PxProng2, a3_hf_cand::PyProng2>,                                       \
    a3_hf_cand::PVectorProng2<a3_hf_cand::PxProng2, a3_hf_cand::PyProng2, a3_hf_cand::PzProng2>

namespace a3DecayMap
{
DECLARE_SOA_COLUMN(DecayMap, decayMap, uint32_t); //! simple map to process passing / not passing criteria
} // namespace a3DecayMap
DECLARE_SOA_TABLE(Alice3DecayMaps, "AOD", "ALICE3DECAYMAPS",
                  a3DecayMap::DecayMap);

using Alice3DecayMap = Alice3DecayMaps::iterator;

namespace a3D0meson
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float);  //! positive track
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float);  //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float);  //!
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float);  //! negative track
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float);  //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float);  //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0,  //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::a3D0meson::pxProng0 + 1.f * aod::a3D0meson::pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::a3D0meson::pyProng0 + 1.f * aod::a3D0meson::pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::a3D0meson::pzProng0 + 1.f * aod::a3D0meson::pzProng1);
DECLARE_SOA_COLUMN(Pt, pt, float); //!
DECLARE_SOA_COLUMN(M, m, float);   //!
DECLARE_SOA_DYNAMIC_COLUMN(E, e,   //!
                           [](float px, float py, float pz, double m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_COLUMN(Eta, eta, float); //!
DECLARE_SOA_COLUMN(Phi, phi, float); //!
DECLARE_SOA_COLUMN(Y, y, float);
} // namespace a3D0meson
DECLARE_SOA_TABLE(Alice3D0Meson, "AOD", "ALICE3D0MESON", //!
                  o2::soa::Index<>,
                  a3D0meson::CollisionId,
                  a3D0meson::PxProng0, a3D0meson::PyProng0, a3D0meson::PzProng0, // positive track
                  a3D0meson::PxProng1, a3D0meson::PyProng1, a3D0meson::PzProng1, // negative track
                  a3D0meson::PtProng0<a3D0meson::PxProng0, a3D0meson::PyProng0>,
                  a3D0meson::PtProng1<a3D0meson::PxProng1, a3D0meson::PyProng1>,
                  a3D0meson::Px, a3D0meson::Py, a3D0meson::Pz,
                  a3D0meson::Pt,
                  a3D0meson::M,
                  a3D0meson::E<a3D0meson::Px, a3D0meson::Py, a3D0meson::Pz, a3D0meson::M>,
                  a3D0meson::Eta,
                  a3D0meson::Phi,
                  a3D0meson::Y);

namespace a3D0Selection
{
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);       //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int); //!
} // namespace a3D0Selection
DECLARE_SOA_TABLE(Alice3D0Sel, "AOD", "ALICE3D0SEL", //!
                  a3D0Selection::IsSelD0,
                  a3D0Selection::IsSelD0bar);

namespace a3D0MCTruth
{
DECLARE_SOA_COLUMN(McTruthInfo, mcTruthInfo, int); //! 0 for bkg, 1 for true D0, 2 for true D0bar
} // namespace a3D0MCTruth
DECLARE_SOA_TABLE(Alice3D0MCTruth, "AOD", "ALICE3D0MCTRUTH", //!
                  a3D0MCTruth::McTruthInfo);                 //!

// Now let's define the Lc to pKpi table
namespace a3_hf_cand_3prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
DECLARE_SOA_COLUMN(Px, px, float);              //!
DECLARE_SOA_COLUMN(Py, py, float);              //!
DECLARE_SOA_COLUMN(Pz, pz, float);              //!
DECLARE_SOA_COLUMN(Pt, pt, float);              //!
DECLARE_SOA_DYNAMIC_COLUMN(M, m,
                           [](float px0, float py0, float pz0,
                              float px1, float py1, float pz1,
                              float px2, float py2, float pz2,
                              const std::array<double, 3>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0},
                                                                                                        std::array{px1, py1, pz1},
                                                                                                        std::array{px2, py2, pz2}},
                                                                                             m); });
DECLARE_SOA_DYNAMIC_COLUMN(E, e, //!
                           [](float px, float py, float pz, const float m) -> float { return RecoDecay::e(px, py, pz, m); });
DECLARE_SOA_COLUMN(Eta, eta, float); //!
DECLARE_SOA_COLUMN(Phi, phi, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(Y, y,
                           [](float px, float py, float pz, const float m) -> float { return RecoDecay::y(std::array{px, py, pz}, m); });
} // namespace a3_hf_cand_3prong
DECLARE_SOA_TABLE(Alice3Cand3Ps, "AOD", "ALICE3CAND3P", //!
                  o2::soa::Index<>,
                  // general candidate properties
                  HFCAND_COLUMNS,
                  HFPRONG0_COLUMNS,
                  HFPRONG1_COLUMNS,
                  HFPRONG2_COLUMNS,
                  // candidate kinematics
                  a3_hf_cand_3prong::Eta,
                  a3_hf_cand_3prong::Phi,
                  a3_hf_cand_3prong::Pt,
                  // prong properties
                  a3_hf_cand::PxProng0, a3_hf_cand::PyProng0, a3_hf_cand::PzProng0, // proton track
                  a3_hf_cand::PxProng1, a3_hf_cand::PyProng1, a3_hf_cand::PzProng1, // kaon track
                  a3_hf_cand::PxProng2, a3_hf_cand::PyProng2, a3_hf_cand::PzProng2, // pion track
                  a3_hf_cand::ImpactParameterY0, a3_hf_cand::ImpactParameterY1, a3_hf_cand::ImpactParameterY2,
                  a3_hf_cand::ErrorImpactParameterY0, a3_hf_cand::ErrorImpactParameterY1, a3_hf_cand::ErrorImpactParameterY2,
                  a3_hf_cand::ImpactParameterZ0, a3_hf_cand::ImpactParameterZ1, a3_hf_cand::ImpactParameterZ2,
                  a3_hf_cand::ErrorImpactParameterZ0, a3_hf_cand::ErrorImpactParameterZ1, a3_hf_cand::ErrorImpactParameterZ2,
                  // Candidate momenta
                  a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py, a3_hf_cand_3prong::Pz,
                  // dynamic candidate properties
                  a3_hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ZSecondaryVertex, a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py, a3_hf_cand_3prong::Pz>,
                  a3_hf_cand::CpaXY<collision::PosX, collision::PosY, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py>,
                  a3_hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, a3_hf_cand::XSecondaryVertex, a3_hf_cand::YSecondaryVertex, a3_hf_cand::ZSecondaryVertex, a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py, a3_hf_cand_3prong::Pz>,
                  a3_hf_cand_3prong::Y<a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py, a3_hf_cand_3prong::Pz>,
                  a3_hf_cand_3prong::M<a3_hf_cand::PxProng0, a3_hf_cand::PyProng0, a3_hf_cand::PzProng0,
                                       a3_hf_cand::PxProng1, a3_hf_cand::PyProng1, a3_hf_cand::PzProng1,
                                       a3_hf_cand::PxProng2, a3_hf_cand::PyProng2, a3_hf_cand::PzProng2>,
                  a3_hf_cand_3prong::E<a3_hf_cand_3prong::Px, a3_hf_cand_3prong::Py, a3_hf_cand_3prong::Pz>);

namespace a3_hf_sel_3prong
{
DECLARE_SOA_COLUMN(IsSelMassHypo0, isSelMassHypo0, bool); //!
DECLARE_SOA_COLUMN(IsSelMassHypo1, isSelMassHypo1, bool); //!

// PID selection
enum PidSels {
  None = 0,
  TrkProng0,
  RichProng0,
  InnTofProng0,
  OutTofProng0,
  TrkProng1,
  RichProng1,
  InnTofProng1,
  OutTofProng1,
  TrkProng2,
  RichProng2,
  InnTofProng2,
  OutTofProng2,
  NPidSelections
};
DECLARE_SOA_COLUMN(PidBitMask, pidBitMask, uint32_t); //!
} // namespace a3_hf_sel_3prong
DECLARE_SOA_TABLE(Alice3Sel3Ps, "AOD", "ALICE3SEL3P", //!
                  a3_hf_sel_3prong::IsSelMassHypo0,
                  a3_hf_sel_3prong::IsSelMassHypo1,
                  a3_hf_sel_3prong::PidBitMask);

namespace a3_mc_truth
{
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int);           //!
DECLARE_SOA_COLUMN(BHadMotherPtRec, bHadMotherPtRec, float); //!
DECLARE_SOA_COLUMN(FlagMcRec, flagMcRec, int);               //!
DECLARE_SOA_COLUMN(ParticleMcRec, particleMcRec, int);       //!
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int);           //!
DECLARE_SOA_COLUMN(BHadMotherPtGen, bHadMotherPtGen, float); //!
DECLARE_SOA_COLUMN(FlagMcGen, flagMcGen, int);               //!
} // namespace a3_mc_truth
DECLARE_SOA_TABLE(Alice3McRecFlags, "AOD", "ALICE3MCRECFLAG", //!
                  a3_mc_truth::OriginMcRec,
                  a3_mc_truth::BHadMotherPtRec,
                  a3_mc_truth::FlagMcRec,
                  a3_mc_truth::ParticleMcRec);

DECLARE_SOA_TABLE(Alice3McGenFlags, "AOD", "ALICE3MCGENFLAG", //!
                  a3_mc_truth::OriginMcGen,
                  a3_mc_truth::BHadMotherPtGen,
                  a3_mc_truth::FlagMcGen);

DECLARE_SOA_TABLE(Alice3PidLcs, "AOD", "ALICE3PIDLC", //!
                  a3_hf_cand::NSigTrkPr0,
                  a3_hf_cand::NSigRichPr0,
                  a3_hf_cand::NSigInnTofPr0,
                  a3_hf_cand::NSigOutTofPr0,
                  a3_hf_cand::NSigTrkKa1,
                  a3_hf_cand::NSigRichKa1,
                  a3_hf_cand::NSigInnTofKa1,
                  a3_hf_cand::NSigOutTofKa1,
                  a3_hf_cand::NSigTrkPi2,
                  a3_hf_cand::NSigRichPi2,
                  a3_hf_cand::NSigInnTofPi2,
                  a3_hf_cand::NSigOutTofPi2);

DECLARE_SOA_TABLE(Alice3Ml3Ps, "AOD", "ALICE3ML3P", //!
                  a3_hf_cand::MlScore0,
                  a3_hf_cand::MlScore1,
                  a3_hf_cand::MlScore2);

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_
