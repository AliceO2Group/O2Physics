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

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otfcascade
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                     //!
DECLARE_SOA_INDEX_COLUMN_FULL(CascadeTrack, cascadeTrack, int, Tracks, "_Cascade"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(BachTrack, bachTrack, int, Tracks, "_Bach");          //!

// topo vars
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DCACascadeDaughters, dcaCascadeDaughters, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
DECLARE_SOA_COLUMN(CascRadiusMC, cascRadiusMC, float);
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(MXi, mXi, float);

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
                  otfcascade::BachTrackId,
                  otfcascade::DCAV0Daughters,
                  otfcascade::DCACascadeDaughters,
                  otfcascade::V0Radius,
                  otfcascade::CascRadius,
                  otfcascade::CascRadiusMC,
                  otfcascade::MLambda,
                  otfcascade::MXi,
                  otfcascade::FindableClusters,
                  otfcascade::FoundClusters);

using UpgradeCascade = UpgradeCascades::iterator;

namespace otfv0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

// topo vars
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0Daughters, float);
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
                  otfv0::PosTrackId,
                  otfv0::NegTrackId,
                  otfv0::DCAV0Daughters,
                  otfv0::V0Radius,
                  otfv0::MLambda,
                  otfv0::MAntiLambda,
                  otfv0::MK0,
                  otfv0::Pt);

using UpgradeV0 = UpgradeV0s::iterator;
} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFSTRANGENESS_H_
