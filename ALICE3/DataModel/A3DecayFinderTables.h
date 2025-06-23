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
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"

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

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_
