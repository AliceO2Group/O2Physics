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
#ifndef PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_
#define PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

//===========================================================================
// For aiding in building: tag those candidates that are interesting
namespace o2::aod
{
namespace hyhe4tag
{
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?
// MC association bools
DECLARE_SOA_COLUMN(IsTrueHyHe4, isTrueHyHe4, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiHyHe4, isTrueAntiHyHe4, bool); //! PDG checked correctly in MC
// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsHyHe4Candidate, isHyHe4Candidate, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsAntiHyHe4Candidate, isAntiHyHe4Candidate, bool); //! compatible with dE/dx hypotheses
} // namespace hyhe4tag
DECLARE_SOA_TABLE(HyHe4Tags, "AOD", "HYHE4TAGS",
                  hyhe4tag::IsInteresting,
                  hyhe4tag::IsTrueHyHe4,
                  hyhe4tag::IsTrueAntiHyHe4,
                  hyhe4tag::IsHyHe4Candidate,
                  hyhe4tag::IsAntiHyHe4Candidate);

//===========================================================================
// The actual analysis data model for hyperhelium4
namespace hyhe4data
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(Helium3Track, prong0Track, int, Tracks, "_Helium3"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(ProtonTrack, prong1Track, int, Tracks, "_Proton");   //!
DECLARE_SOA_INDEX_COLUMN_FULL(PionTrack, prong2Track, int, Tracks, "_Pion");       //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                    //!
DECLARE_SOA_INDEX_COLUMN(Decay3Body, decay3Bodys);                                 //!

// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(Sign, sign, int);             //! sign (positive = matter)
DECLARE_SOA_COLUMN(Helium3X, helium3X, float);   //! prong0 track X at min
DECLARE_SOA_COLUMN(ProtonX, protonX, float);     //! prong1 track X at min
DECLARE_SOA_COLUMN(PionX, pionX, float);         //! prong2 track X at min
DECLARE_SOA_COLUMN(PxHelium3, pxHelium3, float); //! prong 0 track px at min; NB already charge-corrected! Not rigidity anymore!
DECLARE_SOA_COLUMN(PyHelium3, pyHelium3, float); //! prong 0 track py at min; NB already charge-corrected! Not rigidity anymore!
DECLARE_SOA_COLUMN(PzHelium3, pzHelium3, float); //! prong 0 track pz at min; NB already charge-corrected! Not rigidity anymore!
DECLARE_SOA_COLUMN(PxProton, pxProton, float);   //! prong 1 track px at min
DECLARE_SOA_COLUMN(PyProton, pyProton, float);   //! prong 1 track py at min
DECLARE_SOA_COLUMN(PzProton, pzProton, float);   //! prong 1 track pz at min
DECLARE_SOA_COLUMN(PxPion, pxPion, float);       //! prong 1 track px at min
DECLARE_SOA_COLUMN(PyPion, pyPion, float);       //! prong 1 track py at min
DECLARE_SOA_COLUMN(PzPion, pzPion, float);       //! prong 1 track pz at min
DECLARE_SOA_COLUMN(X, x, float);                 //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);                 //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);                 //! decay position Z

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCADaughters, dcaDaughters, float);     //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCAHelium3ToPV, dcaHelium3ToPV, float); //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCAProtonToPV, dcaProtonToPV, float);   //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCAPionToPV, dcaPionToPV, float);       //! DCA positive prong to PV

// DCA to primary vertex
DECLARE_SOA_COLUMN(DCAxyHyHe4ToPV, dcaxyHyHe4ToPV, float); //! DCAxy
DECLARE_SOA_COLUMN(DCAzHyHe4ToPV, dcazHyHe4ToPV, float);   //! DCAz

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! pT
                           [](float pxp0, float pyp0, float pxp1, float pyp1, float pxp2, float pyp2) -> float { return RecoDecay::sqrtSumOfSquares(pxp0 + pxp1 + pxp2, pyp0 + pyp1 + pyp2); });
DECLARE_SOA_DYNAMIC_COLUMN(PtHelium3, ptHelium3, //! pT of prong 0 (identified as Helium-3)
                           [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtProton, ptProton, //! pT of prong 0 (identified as proton)
                           [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(PtPion, ptPion, //! pT of prong 0 (identified as proton)
                           [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });

DECLARE_SOA_DYNAMIC_COLUMN(DecayRadius, decayRadius, //! decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! mass under hyperhelium-4 hypo
                           [](float pxp0, float pyp0, float pzp0, float pxp1, float pyp1, float pzp1, float pxp2, float pyp2, float pzp2) -> float { return RecoDecay::m(std::array{std::array{pxp0, pyp0, pzp0}, std::array{pxp1, pyp1, pzp1}, std::array{pxp2, pyp2, pzp2}}, std::array{o2::constants::physics::MassHelium3, o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(YHyHe4, yHyHe4, //! y -> FIXME add Hyperhelium4 mass to physics constants
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, 3.929); });

// Standard expression columns - note correct momentum for the He3 daughter
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! px
                              float, 1.f * aod::hyhe4data::pxHelium3 + 1.f * aod::hyhe4data::pxProton + 1.f * aod::hyhe4data::pxPion);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! py
                              float, 1.f * aod::hyhe4data::pyHelium3 + 1.f * aod::hyhe4data::pyProton + 1.f * aod::hyhe4data::pyPion);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! pz
                              float, 1.f * aod::hyhe4data::pzHelium3 + 1.f * aod::hyhe4data::pzProton + 1.f * aod::hyhe4data::pzPion);
} // namespace hyhe4data
DECLARE_SOA_TABLE_FULL(StoredHyHe4Datas, "HyHe4Datas", "AOD", "HYHE4DATA", //!
                       o2::soa::Index<>, hyhe4data::Helium3TrackId, hyhe4data::ProtonTrackId, hyhe4data::PionTrackId,
                       hyhe4data::CollisionId, hyhe4data::Decay3BodyId, hyhe4data::Sign,
                       hyhe4data::Helium3X, hyhe4data::ProtonX, hyhe4data::PionX,
                       hyhe4data::PxHelium3, hyhe4data::PyHelium3, hyhe4data::PzHelium3,
                       hyhe4data::PxProton, hyhe4data::PyProton, hyhe4data::PzProton,
                       hyhe4data::PxPion, hyhe4data::PyPion, hyhe4data::PzPion,
                       hyhe4data::X, hyhe4data::Y, hyhe4data::Z,
                       hyhe4data::DCADaughters, hyhe4data::DCAHelium3ToPV, hyhe4data::DCAProtonToPV, hyhe4data::DCAPionToPV,
                       hyhe4data::DCAxyHyHe4ToPV, hyhe4data::DCAzHyHe4ToPV,

                       // dynamic columns
                       hyhe4data::Pt<hyhe4data::PxHelium3, hyhe4data::PyHelium3, hyhe4data::PxProton, hyhe4data::PyProton, hyhe4data::PxPion, hyhe4data::PyPion>,
                       hyhe4data::PtHelium3<hyhe4data::PxHelium3, hyhe4data::PyHelium3>,
                       hyhe4data::PtProton<hyhe4data::PxProton, hyhe4data::PyProton>,
                       hyhe4data::PtPion<hyhe4data::PxPion, hyhe4data::PyPion>,
                       hyhe4data::DecayRadius<hyhe4data::X, hyhe4data::Y>,
                       hyhe4data::M<hyhe4data::PxHelium3, hyhe4data::PyHelium3, hyhe4data::PzHelium3, hyhe4data::PxProton, hyhe4data::PyProton, hyhe4data::PzProton, hyhe4data::PxPion, hyhe4data::PyPion, hyhe4data::PzPion>,
                       hyhe4data::YHyHe4<hyhe4data::Px, hyhe4data::Py, hyhe4data::Pz>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HyHe4Datas, StoredHyHe4Datas, "HYHE4DATAEXT", //!
                                hyhe4data::Px, hyhe4data::Py, hyhe4data::Pz); // the table name has here to be the one with EXT which is not nice and under study

// iterator
using HyHe4Data = HyHe4Datas::iterator;

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_LFHYHEFOURTABLES_H_
