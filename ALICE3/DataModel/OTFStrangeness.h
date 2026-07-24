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
/// \file   OTFStrangeness.h
/// \author David Dobrigkeit Chinellato
/// \since  05/08/2024
/// \brief  Set of tables for the ALICE3 strangeness information
///

#ifndef ALICE3_DATAMODEL_OTFSTRANGENESS_H_
#define ALICE3_DATAMODEL_OTFSTRANGENESS_H_

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cmath>

namespace o2::upgrade::pid
{
static constexpr float NoPidSignal = -999.f;
static constexpr float NoPidSignalThreshold = -990.f;
} // namespace o2::upgrade::pid

namespace o2::aod
{
namespace otfcascade
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                     //!
DECLARE_SOA_INDEX_COLUMN_FULL(CascadeTrack, cascadeTrack, int, Tracks, "_Cascade"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, Tracks, "_Bach");            //!

DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaCascadeDaughters, dcaCascadeDaughters, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
DECLARE_SOA_COLUMN(CascRadiusMC, cascRadiusMC, float);
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(MXi, mXi, float);
DECLARE_SOA_COLUMN(DcaXYCascToPV, dcaXYCascToPV, float);
DECLARE_SOA_COLUMN(DcaZCascToPV, dcaZCascToPV, float);

// strangeness tracking
DECLARE_SOA_COLUMN(FindableClusters, findableClusters, int);
DECLARE_SOA_COLUMN(FoundClusters, foundClusters, int);

} // namespace otfcascade
DECLARE_SOA_TABLE(UpgradeCascades, "AOD", "UPGRADECASCADES",
                  o2::soa::Index<>,
                  otfcascade::CollisionId,
                  otfcascade::CascadeTrackId,
                  otfcascade::PosTrackId,
                  otfcascade::NegTrackId,
                  otfcascade::BachelorId,
                  otfcascade::DcaV0Daughters,
                  otfcascade::DcaCascadeDaughters,
                  otfcascade::V0Radius,
                  otfcascade::CascRadius,
                  otfcascade::CascRadiusMC,
                  otfcascade::MLambda,
                  otfcascade::MXi,
                  otfcascade::FindableClusters,
                  otfcascade::FoundClusters,
                  otfcascade::DcaXYCascToPV,
                  otfcascade::DcaZCascToPV);

using UpgradeCascade = UpgradeCascades::iterator;

DECLARE_SOA_TABLE(A3CascadeMcLabels, "AOD", "A3CASCADEMCLABELS",
                  o2::soa::Index<>, otfcascade::McParticleId);
DECLARE_SOA_TABLE(UpgradeCascadeMcLabels, "AOD", "UPGRADECASCMCLAB",
                  o2::soa::Index<>, otfcascade::McParticleId);

namespace casc_pid
{

// Xi
DECLARE_SOA_COLUMN(NSigmaInnerTofXiBachPi, nSigmaInnerTofXiBachPi, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofXiPosPr, nSigmaInnerTofXiPosPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofXiPosPi, nSigmaInnerTofXiPosPi, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofXiNegPr, nSigmaInnerTofXiNegPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofXiNegPi, nSigmaInnerTofXiNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofXiBachPi, hasInnerTofXiBachPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofXiPosPr, hasInnerTofXiPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofXiPosPi, hasInnerTofXiPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofXiNegPr, hasInnerTofXiNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofXiNegPi, hasInnerTofXiNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(NSigmaOuterTofXiBachPi, nSigmaOuterTofXiBachPi, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofXiPosPr, nSigmaOuterTofXiPosPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofXiPosPi, nSigmaOuterTofXiPosPi, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofXiNegPr, nSigmaOuterTofXiNegPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofXiNegPi, nSigmaOuterTofXiNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofXiBachPi, hasOuterTofXiBachPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofXiPosPr, hasOuterTofXiPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofXiPosPi, hasOuterTofXiPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofXiNegPr, hasOuterTofXiNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofXiNegPi, hasOuterTofXiNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(ExpectedInnerTofXiBachPi, expectedInnerTofXiBachPi, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofXiPosPr, expectedInnerTofXiPosPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofXiPosPi, expectedInnerTofXiPosPi, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofXiNegPr, expectedInnerTofXiNegPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofXiNegPi, expectedInnerTofXiNegPi, float);

DECLARE_SOA_COLUMN(ExpectedOuterTofXiBachPi, expectedOuterTofXiBachPi, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofXiPosPr, expectedOuterTofXiPosPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofXiPosPi, expectedOuterTofXiPosPi, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofXiNegPr, expectedOuterTofXiNegPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofXiNegPi, expectedOuterTofXiNegPi, float);

DECLARE_SOA_COLUMN(MeasuredInnerTofXiBachPi, measuredInnerTofXiBachPi, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofXiPosPr, measuredInnerTofXiPosPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofXiPosPi, measuredInnerTofXiPosPi, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofXiNegPr, measuredInnerTofXiNegPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofXiNegPi, measuredInnerTofXiNegPi, float);

DECLARE_SOA_COLUMN(MeasuredOuterTofXiBachPi, measuredOuterTofXiBachPi, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofXiPosPr, measuredOuterTofXiPosPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofXiPosPi, measuredOuterTofXiPosPi, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofXiNegPr, measuredOuterTofXiNegPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofXiNegPi, measuredOuterTofXiNegPi, float);

// Omega
DECLARE_SOA_COLUMN(NSigmaInnerTofOmegaBachKa, nSigmaInnerTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofOmegaPosPr, nSigmaInnerTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofOmegaPosPi, nSigmaInnerTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofOmegaNegPr, nSigmaInnerTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofOmegaNegPi, nSigmaInnerTofOmegaNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofOmegaBachKa, hasInnerTofOmegaBachKa,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofOmegaPosPr, hasInnerTofOmegaPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofOmegaPosPi, hasInnerTofOmegaPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofOmegaNegPr, hasInnerTofOmegaNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofOmegaNegPi, hasInnerTofOmegaNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(NSigmaOuterTofOmegaBachKa, nSigmaOuterTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofOmegaPosPr, nSigmaOuterTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofOmegaPosPi, nSigmaOuterTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofOmegaNegPr, nSigmaOuterTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofOmegaNegPi, nSigmaOuterTofOmegaNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofOmegaBachKa, hasOuterTofOmegaBachKa,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofOmegaPosPr, hasOuterTofOmegaPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofOmegaPosPi, hasOuterTofOmegaPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofOmegaNegPr, hasOuterTofOmegaNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofOmegaNegPi, hasOuterTofOmegaNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(ExpectedInnerTofOmegaBachKa, expectedInnerTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofOmegaPosPr, expectedInnerTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofOmegaPosPi, expectedInnerTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofOmegaNegPr, expectedInnerTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofOmegaNegPi, expectedInnerTofOmegaNegPi, float);

DECLARE_SOA_COLUMN(ExpectedOuterTofOmegaBachKa, expectedOuterTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofOmegaPosPr, expectedOuterTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofOmegaPosPi, expectedOuterTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofOmegaNegPr, expectedOuterTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofOmegaNegPi, expectedOuterTofOmegaNegPi, float);

DECLARE_SOA_COLUMN(MeasuredInnerTofOmegaBachKa, measuredInnerTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofOmegaPosPr, measuredInnerTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofOmegaPosPi, measuredInnerTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofOmegaNegPr, measuredInnerTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofOmegaNegPi, measuredInnerTofOmegaNegPi, float);

DECLARE_SOA_COLUMN(MeasuredOuterTofOmegaBachKa, measuredOuterTofOmegaBachKa, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofOmegaPosPr, measuredOuterTofOmegaPosPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofOmegaPosPi, measuredOuterTofOmegaPosPi, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofOmegaNegPr, measuredOuterTofOmegaNegPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofOmegaNegPi, measuredOuterTofOmegaNegPi, float);

} // namespace casc_pid

DECLARE_SOA_TABLE(A3XiInnerTofPid, "AOD", "A3XIITOFPID",
                  casc_pid::NSigmaInnerTofXiBachPi,
                  casc_pid::NSigmaInnerTofXiPosPr,
                  casc_pid::NSigmaInnerTofXiPosPi,
                  casc_pid::NSigmaInnerTofXiNegPr,
                  casc_pid::NSigmaInnerTofXiNegPi,
                  casc_pid::HasInnerTofXiBachPi<casc_pid::NSigmaInnerTofXiBachPi>,
                  casc_pid::HasInnerTofXiPosPr<casc_pid::NSigmaInnerTofXiPosPr>,
                  casc_pid::HasInnerTofXiPosPi<casc_pid::NSigmaInnerTofXiPosPi>,
                  casc_pid::HasInnerTofXiNegPr<casc_pid::NSigmaInnerTofXiNegPr>,
                  casc_pid::HasInnerTofXiNegPi<casc_pid::NSigmaInnerTofXiNegPi>);

DECLARE_SOA_TABLE(A3XiOuterTofPid, "AOD", "A3XIOTOFPID",
                  casc_pid::NSigmaOuterTofXiBachPi,
                  casc_pid::NSigmaOuterTofXiPosPr,
                  casc_pid::NSigmaOuterTofXiPosPi,
                  casc_pid::NSigmaOuterTofXiNegPr,
                  casc_pid::NSigmaOuterTofXiNegPi,
                  casc_pid::HasOuterTofXiBachPi<casc_pid::NSigmaOuterTofXiBachPi>,
                  casc_pid::HasOuterTofXiPosPr<casc_pid::NSigmaOuterTofXiPosPr>,
                  casc_pid::HasOuterTofXiPosPi<casc_pid::NSigmaOuterTofXiPosPi>,
                  casc_pid::HasOuterTofXiNegPr<casc_pid::NSigmaOuterTofXiNegPr>,
                  casc_pid::HasOuterTofXiNegPi<casc_pid::NSigmaOuterTofXiNegPi>);

DECLARE_SOA_TABLE(A3OmegaInnerTofPid, "AOD", "A3OMITOFPID",
                  casc_pid::NSigmaInnerTofOmegaBachKa,
                  casc_pid::NSigmaInnerTofOmegaPosPr,
                  casc_pid::NSigmaInnerTofOmegaPosPi,
                  casc_pid::NSigmaInnerTofOmegaNegPr,
                  casc_pid::NSigmaInnerTofOmegaNegPi,
                  casc_pid::HasInnerTofOmegaBachKa<casc_pid::NSigmaInnerTofOmegaBachKa>,
                  casc_pid::HasInnerTofOmegaPosPr<casc_pid::NSigmaInnerTofOmegaPosPr>,
                  casc_pid::HasInnerTofOmegaPosPi<casc_pid::NSigmaInnerTofOmegaPosPi>,
                  casc_pid::HasInnerTofOmegaNegPr<casc_pid::NSigmaInnerTofOmegaNegPr>,
                  casc_pid::HasInnerTofOmegaNegPi<casc_pid::NSigmaInnerTofOmegaNegPi>);

DECLARE_SOA_TABLE(A3OmegaOuterTofPid, "AOD", "A3OMOTOFPID",
                  casc_pid::NSigmaOuterTofOmegaBachKa,
                  casc_pid::NSigmaOuterTofOmegaPosPr,
                  casc_pid::NSigmaOuterTofOmegaPosPi,
                  casc_pid::NSigmaOuterTofOmegaNegPr,
                  casc_pid::NSigmaOuterTofOmegaNegPi,
                  casc_pid::HasOuterTofOmegaBachKa<casc_pid::NSigmaOuterTofOmegaBachKa>,
                  casc_pid::HasOuterTofOmegaPosPr<casc_pid::NSigmaOuterTofOmegaPosPr>,
                  casc_pid::HasOuterTofOmegaPosPi<casc_pid::NSigmaOuterTofOmegaPosPi>,
                  casc_pid::HasOuterTofOmegaNegPr<casc_pid::NSigmaOuterTofOmegaNegPr>,
                  casc_pid::HasOuterTofOmegaNegPi<casc_pid::NSigmaOuterTofOmegaNegPi>);

DECLARE_SOA_TABLE(A3XiExpectedInnerTimes, "AOD", "A3XIITIMES",
                  casc_pid::ExpectedInnerTofXiBachPi,
                  casc_pid::ExpectedInnerTofXiPosPr,
                  casc_pid::ExpectedInnerTofXiPosPi,
                  casc_pid::ExpectedInnerTofXiNegPr,
                  casc_pid::ExpectedInnerTofXiNegPi,
                  casc_pid::MeasuredInnerTofXiBachPi,
                  casc_pid::MeasuredInnerTofXiPosPr,
                  casc_pid::MeasuredInnerTofXiPosPi,
                  casc_pid::MeasuredInnerTofXiNegPr,
                  casc_pid::MeasuredInnerTofXiNegPi);

DECLARE_SOA_TABLE(A3XiExpectedOuterTimes, "AOD", "A3XIOTIMES",
                  casc_pid::ExpectedOuterTofXiBachPi,
                  casc_pid::ExpectedOuterTofXiPosPr,
                  casc_pid::ExpectedOuterTofXiPosPi,
                  casc_pid::ExpectedOuterTofXiNegPr,
                  casc_pid::ExpectedOuterTofXiNegPi,
                  casc_pid::MeasuredOuterTofXiBachPi,
                  casc_pid::MeasuredOuterTofXiPosPr,
                  casc_pid::MeasuredOuterTofXiPosPi,
                  casc_pid::MeasuredOuterTofXiNegPr,
                  casc_pid::MeasuredOuterTofXiNegPi);

DECLARE_SOA_TABLE(A3OmegaExpectedInnerTimes, "AOD", "A3OMITIMES",
                  casc_pid::ExpectedInnerTofOmegaBachKa,
                  casc_pid::ExpectedInnerTofOmegaPosPr,
                  casc_pid::ExpectedInnerTofOmegaPosPi,
                  casc_pid::ExpectedInnerTofOmegaNegPr,
                  casc_pid::ExpectedInnerTofOmegaNegPi,
                  casc_pid::MeasuredInnerTofOmegaBachKa,
                  casc_pid::MeasuredInnerTofOmegaPosPr,
                  casc_pid::MeasuredInnerTofOmegaPosPi,
                  casc_pid::MeasuredInnerTofOmegaNegPr,
                  casc_pid::MeasuredInnerTofOmegaNegPi);

DECLARE_SOA_TABLE(A3OmegaExpectedOuterTimes, "AOD", "A3OMOTIMES",
                  casc_pid::ExpectedOuterTofOmegaBachKa,
                  casc_pid::ExpectedOuterTofOmegaPosPr,
                  casc_pid::ExpectedOuterTofOmegaPosPi,
                  casc_pid::ExpectedOuterTofOmegaNegPr,
                  casc_pid::ExpectedOuterTofOmegaNegPi,
                  casc_pid::MeasuredOuterTofOmegaBachKa,
                  casc_pid::MeasuredOuterTofOmegaPosPr,
                  casc_pid::MeasuredOuterTofOmegaPosPi,
                  casc_pid::MeasuredOuterTofOmegaNegPr,
                  casc_pid::MeasuredOuterTofOmegaNegPi);

namespace otfv0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //! index of the mc particle corresponding to the V0

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);
DECLARE_SOA_COLUMN(MK0, mK0, float);

// kinematics
DECLARE_SOA_COLUMN(Pt, pt, float);

} // namespace otfv0
DECLARE_SOA_TABLE(UpgradeV0s, "AOD", "UPGRADEV0S",
                  o2::soa::Index<>,
                  otfv0::CollisionId,
                  otfv0::V0Id,
                  otfv0::PosTrackId,
                  otfv0::NegTrackId,
                  otfv0::DcaV0Daughters,
                  otfv0::V0Radius,
                  otfv0::MLambda,
                  otfv0::MAntiLambda,
                  otfv0::MK0,
                  otfv0::Pt);

using UpgradeV0 = UpgradeV0s::iterator;

namespace candidatev0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

// Label to MC particle
DECLARE_SOA_INDEX_COLUMN_FULL(McParticle, mcParticle, int, McParticles, ""); //! label to the MC particle corresponding to the V0

// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(PxPos, pxPos, float); //! positive track px at min
DECLARE_SOA_COLUMN(PyPos, pyPos, float); //! positive track py at min
DECLARE_SOA_COLUMN(PzPos, pzPos, float); //! positive track pz at min
DECLARE_SOA_COLUMN(PxNeg, pxNeg, float); //! negative track px at min
DECLARE_SOA_COLUMN(PyNeg, pyNeg, float); //! negative track py at min
DECLARE_SOA_COLUMN(PzNeg, pzNeg, float); //! negative track pz at min
DECLARE_SOA_COLUMN(X, x, float);         //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);         //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);         //! decay position Z

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);

//______________________________________________________
// DYNAMIC COLUMNS

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! V0 px
                           [](float pxPos, float pxNeg) -> float { return pxPos + pxNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! V0 py
                           [](float pyPos, float pyNeg) -> float { return pyPos + pyNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! V0 pz
                           [](float pzPos, float pzNeg) -> float { return pzPos + pzNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Transverse momentum in GeV/c
                           [](float pxPos, float pyPos, float pxNeg, float pyNeg) -> float {
                             return RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Total momentum in GeV/c
                           [](float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg) -> float {
                             return RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Phi in the range [0, 2pi)
                           [](float pxPos, float pyPos, float pxNeg, float pyNeg) -> float { return RecoDecay::phi(pxPos + pxNeg, pyPos + pyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg) -> float {
                             return RecoDecay::eta(std::array{pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg});
                           });
// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0radius, v0radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distOverTotMom, //! PV to V0decay distance over total momentum
                           [](float X, float Y, float Z, float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg, float pvX, float pvY, float pvZ) {
                             float p = RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (p + 1E-10);
                           });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtArm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });
// Mass assumption
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                           });
// Rapidity
DECLARE_SOA_DYNAMIC_COLUMN(YK0Short, yK0Short, //! V0 y with K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::y(std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}, o2::constants::physics::MassKaonNeutral);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(YLambda, yLambda, //! V0 y with lambda or antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::y(std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}, o2::constants::physics::MassLambda);
                           });
// Daughter track momenta
DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativePt, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivePt, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeEta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(std::array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativePhi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveEta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivePhi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });
} // namespace candidatev0
DECLARE_SOA_TABLE(V0CandidateIndices, "AOD", "V0CANDIDATEINDEX", //! index table
                  o2::soa::Index<>, candidatev0::CollisionId, candidatev0::PosTrackId, candidatev0::NegTrackId, candidatev0::McParticleId);

DECLARE_SOA_TABLE(V0CandidateCores, "AOD", "V0CANDIDATECORE",
                  o2::soa::Index<>,
                  candidatev0::X, candidatev0::Y, candidatev0::Z,
                  candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos,
                  candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg,
                  candidatev0::DcaV0Daughters, candidatev0::DcaPosToPV, candidatev0::DcaNegToPV,
                  candidatev0::CosPA, candidatev0::DcaV0ToPV,
                  candidatev0::Px<candidatev0::PxPos, candidatev0::PxNeg>,
                  candidatev0::Py<candidatev0::PyPos, candidatev0::PyNeg>,
                  candidatev0::Pz<candidatev0::PzPos, candidatev0::PzNeg>,
                  candidatev0::Pt<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::P<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::Phi<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::Eta<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::V0radius<candidatev0::X, candidatev0::Y>,
                  candidatev0::DistOverTotMom<candidatev0::X, candidatev0::Y, candidatev0::Z, candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::Alpha<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::QtArm<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MAntiLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MK0Short<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::YK0Short<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::YLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::NegativePt<candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::PositivePt<candidatev0::PxPos, candidatev0::PyPos>,
                  candidatev0::NegativeEta<candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::NegativePhi<candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::PositiveEta<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos>,
                  candidatev0::PositivePhi<candidatev0::PxPos, candidatev0::PyPos>);

using V0CandidateCore = V0CandidateCores::iterator;

namespace v0_pid
{

// K0S
DECLARE_SOA_COLUMN(NSigmaInnerTofK0SPosPi, nSigmaInnerTofK0SPosPi, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofK0SNegPi, nSigmaInnerTofK0SNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofK0SPosPi, hasInnerTofK0SPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofK0SNegPi, hasInnerTofK0SNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(NSigmaOuterTofK0SPosPi, nSigmaOuterTofK0SPosPi, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofK0SNegPi, nSigmaOuterTofK0SNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofK0SPosPi, hasOuterTofK0SPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofK0SNegPi, hasOuterTofK0SNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(ExpectedInnerTofK0SPosPi, expectedInnerTofK0SPosPi, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofK0SNegPi, expectedInnerTofK0SNegPi, float);

DECLARE_SOA_COLUMN(ExpectedOuterTofK0SPosPi, expectedOuterTofK0SPosPi, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofK0SNegPi, expectedOuterTofK0SNegPi, float);

DECLARE_SOA_COLUMN(MeasuredInnerTofK0SPosPi, measuredInnerTofK0SPosPi, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofK0SNegPi, measuredInnerTofK0SNegPi, float);

DECLARE_SOA_COLUMN(MeasuredOuterTofK0SPosPi, measuredOuterTofK0SPosPi, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofK0SNegPi, measuredOuterTofK0SNegPi, float);

// Lambda
DECLARE_SOA_COLUMN(NSigmaInnerTofLambdaPosPr, nSigmaInnerTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofLambdaPosPi, nSigmaInnerTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofLambdaNegPr, nSigmaInnerTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(NSigmaInnerTofLambdaNegPi, nSigmaInnerTofLambdaNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofLambdaPosPr, hasInnerTofLambdaPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofLambdaPosPi, hasInnerTofLambdaPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofLambdaNegPr, hasInnerTofLambdaNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasInnerTofLambdaNegPi, hasInnerTofLambdaNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(NSigmaOuterTofLambdaPosPr, nSigmaOuterTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofLambdaPosPi, nSigmaOuterTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofLambdaNegPr, nSigmaOuterTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(NSigmaOuterTofLambdaNegPi, nSigmaOuterTofLambdaNegPi, float);
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofLambdaPosPr, hasOuterTofLambdaPosPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofLambdaPosPi, hasOuterTofLambdaPosPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofLambdaNegPr, hasOuterTofLambdaNegPr,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });
DECLARE_SOA_DYNAMIC_COLUMN(HasOuterTofLambdaNegPi, hasOuterTofLambdaNegPi,
                           [](float nSigma) -> bool { return (nSigma > o2::upgrade::pid::NoPidSignalThreshold); });

DECLARE_SOA_COLUMN(ExpectedInnerTofLambdaPosPr, expectedInnerTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofLambdaPosPi, expectedInnerTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofLambdaNegPr, expectedInnerTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(ExpectedInnerTofLambdaNegPi, expectedInnerTofLambdaNegPi, float);

DECLARE_SOA_COLUMN(ExpectedOuterTofLambdaPosPr, expectedOuterTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofLambdaPosPi, expectedOuterTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofLambdaNegPr, expectedOuterTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(ExpectedOuterTofLambdaNegPi, expectedOuterTofLambdaNegPi, float);

DECLARE_SOA_COLUMN(MeasuredInnerTofLambdaPosPr, measuredInnerTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofLambdaPosPi, measuredInnerTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofLambdaNegPr, measuredInnerTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(MeasuredInnerTofLambdaNegPi, measuredInnerTofLambdaNegPi, float);

DECLARE_SOA_COLUMN(MeasuredOuterTofLambdaPosPr, measuredOuterTofLambdaPosPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofLambdaPosPi, measuredOuterTofLambdaPosPi, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofLambdaNegPr, measuredOuterTofLambdaNegPr, float);
DECLARE_SOA_COLUMN(MeasuredOuterTofLambdaNegPi, measuredOuterTofLambdaNegPi, float);
} // namespace v0_pid

DECLARE_SOA_TABLE(A3K0SInnerTofPid, "AOD", "A3K0SITOFPID",
                  v0_pid::NSigmaInnerTofK0SPosPi,
                  v0_pid::NSigmaInnerTofK0SNegPi,
                  v0_pid::HasInnerTofK0SPosPi<v0_pid::NSigmaInnerTofK0SPosPi>,
                  v0_pid::HasInnerTofK0SNegPi<v0_pid::NSigmaInnerTofK0SNegPi>);

DECLARE_SOA_TABLE(A3K0SOuterTofPid, "AOD", "A3K0SOTOFPID",
                  v0_pid::NSigmaOuterTofK0SPosPi,
                  v0_pid::NSigmaOuterTofK0SNegPi,
                  v0_pid::HasOuterTofK0SPosPi<v0_pid::NSigmaOuterTofK0SPosPi>,
                  v0_pid::HasOuterTofK0SNegPi<v0_pid::NSigmaOuterTofK0SNegPi>);

DECLARE_SOA_TABLE(A3LambdaInnerTofPid, "AOD", "A3LAITOFPID",
                  v0_pid::NSigmaInnerTofLambdaPosPr,
                  v0_pid::NSigmaInnerTofLambdaPosPi,
                  v0_pid::NSigmaInnerTofLambdaNegPr,
                  v0_pid::NSigmaInnerTofLambdaNegPi,
                  v0_pid::HasInnerTofLambdaPosPr<v0_pid::NSigmaInnerTofLambdaPosPr>,
                  v0_pid::HasInnerTofLambdaPosPi<v0_pid::NSigmaInnerTofLambdaPosPi>,
                  v0_pid::HasInnerTofLambdaNegPr<v0_pid::NSigmaInnerTofLambdaNegPr>,
                  v0_pid::HasInnerTofLambdaNegPi<v0_pid::NSigmaInnerTofLambdaNegPi>);

DECLARE_SOA_TABLE(A3LambdaOuterTofPid, "AOD", "A3LAOTOFPID",
                  v0_pid::NSigmaOuterTofLambdaPosPr,
                  v0_pid::NSigmaOuterTofLambdaPosPi,
                  v0_pid::NSigmaOuterTofLambdaNegPr,
                  v0_pid::NSigmaOuterTofLambdaNegPi,
                  v0_pid::HasOuterTofLambdaPosPr<v0_pid::NSigmaOuterTofLambdaPosPr>,
                  v0_pid::HasOuterTofLambdaPosPi<v0_pid::NSigmaOuterTofLambdaPosPi>,
                  v0_pid::HasOuterTofLambdaNegPr<v0_pid::NSigmaOuterTofLambdaNegPr>,
                  v0_pid::HasOuterTofLambdaNegPi<v0_pid::NSigmaOuterTofLambdaNegPi>);

DECLARE_SOA_TABLE(A3K0SExpectedInnerTimes, "AOD", "A3K0SITIMES",
                  v0_pid::ExpectedInnerTofK0SPosPi,
                  v0_pid::ExpectedInnerTofK0SNegPi,
                  v0_pid::MeasuredInnerTofK0SPosPi,
                  v0_pid::MeasuredInnerTofK0SNegPi);

DECLARE_SOA_TABLE(A3K0SExpectedOuterTimes, "AOD", "A3K0SOTIMES",
                  v0_pid::ExpectedOuterTofK0SPosPi,
                  v0_pid::ExpectedOuterTofK0SNegPi,
                  v0_pid::MeasuredOuterTofK0SPosPi,
                  v0_pid::MeasuredOuterTofK0SNegPi);

DECLARE_SOA_TABLE(A3LambdaExpectedInnerTimes, "AOD", "A3LAITIMES",
                  v0_pid::ExpectedInnerTofLambdaPosPr,
                  v0_pid::ExpectedInnerTofLambdaPosPi,
                  v0_pid::ExpectedInnerTofLambdaNegPr,
                  v0_pid::ExpectedInnerTofLambdaNegPi,
                  v0_pid::MeasuredInnerTofLambdaPosPr,
                  v0_pid::MeasuredInnerTofLambdaPosPi,
                  v0_pid::MeasuredInnerTofLambdaNegPr,
                  v0_pid::MeasuredInnerTofLambdaNegPi);

DECLARE_SOA_TABLE(A3LambdaExpectedOuterTimes, "AOD", "A3LAOTIMES",
                  v0_pid::ExpectedOuterTofLambdaPosPr,
                  v0_pid::ExpectedOuterTofLambdaPosPi,
                  v0_pid::ExpectedOuterTofLambdaNegPr,
                  v0_pid::ExpectedOuterTofLambdaNegPi,
                  v0_pid::MeasuredOuterTofLambdaPosPr,
                  v0_pid::MeasuredOuterTofLambdaPosPi,
                  v0_pid::MeasuredOuterTofLambdaNegPr,
                  v0_pid::MeasuredOuterTofLambdaNegPi);

} // namespace o2::aod
#endif // ALICE3_DATAMODEL_OTFSTRANGENESS_H_
