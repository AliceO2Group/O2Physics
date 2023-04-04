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

/// \brief This task serves to do hadron-(strange hadron) correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : V0s or Cascades
///
/// \author Kai Cui
/// \since
#ifndef O2_ANALYSIS_HSTRANGECORRELATIONTABLES_H_
#define O2_ANALYSIS_HSTRANGECORRELATIONTABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include <cmath>

namespace o2::aod
{
/// _________________________________________
/// Table for storing trigger track indices
namespace triggerTracks
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Trigger"); //!
} // namespace triggerTracks
DECLARE_SOA_TABLE(TriggerTracks, "AOD", "TRIGGERTRACKS", o2::soa::Index<>, triggerTracks::CollisionId, triggerTracks::TrackId);
/// _________________________________________
/// Table for storing associated V0 indices
namespace assocV0s
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_INDEX_COLUMN(V0Data, v0Data);                             //!
DECLARE_SOA_COLUMN(CompatibleK0Short, compatibleK0Short, bool);       // compatible with K0Short
DECLARE_SOA_COLUMN(CompatibleLambda, compatibleLambda, bool);         // compatible with Lambda
DECLARE_SOA_COLUMN(CompatibleAntiLambda, compatibleAntiLambda, bool); // compatible with AntiLambda
} // namespace assocV0s
DECLARE_SOA_TABLE(AssocV0s, "AOD", "ASSOCV0S", o2::soa::Index<>,
                  assocV0s::CollisionId, assocV0s::V0DataId,
                  assocV0s::CompatibleK0Short,
                  assocV0s::CompatibleLambda,
                  assocV0s::CompatibleAntiLambda);
/// _________________________________________
/// Table for storing associated casc indices
namespace assocCascades
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_INDEX_COLUMN(CascData, cascData);                         //!
DECLARE_SOA_COLUMN(CompatibleXiMinus, compatibleXiMinus, bool);       // compatible with XiMinus
DECLARE_SOA_COLUMN(CompatibleXiPlus, compatibleXiPlus, bool);         // compatible with XiPlus
DECLARE_SOA_COLUMN(CompatibleOmegaMinus, compatibleOmegaMinus, bool); // compatible with OmegaMinus
DECLARE_SOA_COLUMN(CompatibleOmegaPlus, compatibleOmegaPlus, bool);   // compatible with OmegaPlus
} // namespace assocCascades
DECLARE_SOA_TABLE(AssocCascades, "AOD", "ASSOCCASCADES", o2::soa::Index<>, assocCascades::CollisionId, assocCascades::CascDataId,
                  assocCascades::CompatibleXiMinus,
                  assocCascades::CompatibleXiPlus,
                  assocCascades::CompatibleOmegaMinus,
                  assocCascades::CompatibleOmegaPlus);
} // namespace o2::aod

#endif // O2_ANALYSIS_HSTRANGECORRELATIONTABLES_H_
