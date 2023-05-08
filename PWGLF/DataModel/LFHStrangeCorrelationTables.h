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
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

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
DECLARE_SOA_COLUMN(MassRegionK0Short, massRegionK0Short, int);        //
DECLARE_SOA_COLUMN(MassRegionLambda, massRegionLambda, int);          //
DECLARE_SOA_COLUMN(MassRegionAntiLambda, massRegionAntiLambda, int);  //
DECLARE_SOA_DYNAMIC_COLUMN(Compatible, compatible,                    //! check compatibility with a hypothesis of a certain number (0 - K0, 1 - L, 2 - Lbar)
                           [](bool cK0Short, bool cLambda, bool cAntiLambda, int value) -> bool {
                             if (value == 0 && cK0Short)
                               return true;
                             if (value == 1 && cLambda)
                               return true;
                             if (value == 2 && cAntiLambda)
                               return true;
                             return false;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(InMassRegionCheck, inMassRegionCheck,
                           [](int rK0Short, int rLambda, int rAntiLambda, int value, int region) -> bool {
                             if (value == 0 && rK0Short == region)
                               return true;
                             if (value == 1 && rLambda == region)
                               return true;
                             if (value == 2 && rAntiLambda == region)
                               return true;
                             return false;
                           });
} // namespace assocV0s
DECLARE_SOA_TABLE(AssocV0s, "AOD", "ASSOCV0S", o2::soa::Index<>,
                  assocV0s::CollisionId, assocV0s::V0DataId,
                  assocV0s::CompatibleK0Short,
                  assocV0s::CompatibleLambda,
                  assocV0s::CompatibleAntiLambda,
                  assocV0s::MassRegionK0Short,
                  assocV0s::MassRegionLambda,
                  assocV0s::MassRegionAntiLambda,
                  assocV0s::Compatible<assocV0s::CompatibleK0Short, assocV0s::CompatibleLambda, assocV0s::CompatibleAntiLambda>,
                  assocV0s::InMassRegionCheck<assocV0s::MassRegionK0Short, assocV0s::MassRegionLambda, assocV0s::MassRegionAntiLambda>);
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
DECLARE_SOA_COLUMN(MassRegionXi, massRegionXi, int);                  //
DECLARE_SOA_COLUMN(MassRegionOmega, massRegionOmega, int);            //
DECLARE_SOA_DYNAMIC_COLUMN(Compatible, compatible,                    //! check compatibility with a hypothesis of a certain number (0 - K0, 1 - L, 2 - Lbar)
                           [](bool cXiMinus, bool cXiPlus, bool cOmegaMinus, bool cOmegaPlus, int value) -> bool {
                             if (value == 0 && cXiMinus)
                               return true;
                             if (value == 1 && cXiPlus)
                               return true;
                             if (value == 2 && cOmegaMinus)
                               return true;
                             if (value == 3 && cOmegaPlus)
                               return true;
                             return false;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(InMassRegionCheck, inMassRegionCheck,
                           [](int rXi, int rOmega, int value, int region) -> bool {
                             if (value == 0 && rXi == region)
                               return true;
                             if (value == 2 && rOmega == region)
                               return true;
                             return false;
                           });
} // namespace assocCascades
DECLARE_SOA_TABLE(AssocCascades, "AOD", "ASSOCCASCADES", o2::soa::Index<>, assocCascades::CollisionId, assocCascades::CascDataId,
                  assocCascades::CompatibleXiMinus,
                  assocCascades::CompatibleXiPlus,
                  assocCascades::CompatibleOmegaMinus,
                  assocCascades::CompatibleOmegaPlus,
                  assocCascades::MassRegionXi,
                  assocCascades::MassRegionOmega,
                  assocCascades::Compatible<assocCascades::CompatibleXiMinus, assocCascades::CompatibleXiPlus, assocCascades::CompatibleOmegaMinus, assocCascades::CompatibleOmegaPlus>,
                  assocCascades::InMassRegionCheck<assocCascades::MassRegionXi, assocCascades::MassRegionOmega>);
} // namespace o2::aod

#endif // O2_ANALYSIS_HSTRANGECORRELATIONTABLES_H_