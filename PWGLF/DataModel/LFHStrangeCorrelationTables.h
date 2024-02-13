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

#ifndef PWGLF_DATAMODEL_LFHSTRANGECORRELATIONTABLES_H_
#define PWGLF_DATAMODEL_LFHSTRANGECORRELATIONTABLES_H_

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
/// Table for storing assoc track indices
namespace assocPions
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                     //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_Assoc"); //!
} // namespace assocPions
DECLARE_SOA_TABLE(AssocPions, "AOD", "ASSOCPIONS", o2::soa::Index<>, assocPions::CollisionId, assocPions::TrackId);
/// _________________________________________
/// Table for storing associated V0 indices
namespace assocV0s
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                       //!
DECLARE_SOA_INDEX_COLUMN(V0Data, v0Data);                             //!
DECLARE_SOA_COLUMN(CompatibleK0Short, compatibleK0Short, bool);       // compatible with K0Short
DECLARE_SOA_COLUMN(CompatibleLambda, compatibleLambda, bool);         // compatible with Lambda
DECLARE_SOA_COLUMN(CompatibleAntiLambda, compatibleAntiLambda, bool); // compatible with AntiLambda
DECLARE_SOA_COLUMN(MCTrueK0Short, mcTrueK0Short, bool);               // true K0Short in MC
DECLARE_SOA_COLUMN(MCTrueLambda, mcTrueLambda, bool);                 // true Lambda in MC
DECLARE_SOA_COLUMN(MCTrueAntiLambda, mcTrueAntiLambda, bool);         // true AntiLambda in MC
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
DECLARE_SOA_DYNAMIC_COLUMN(InvMassRegionCheck, invMassRegionCheck,
                           [](int rK0Short, int rLambda, int rAntiLambda, int value, int region) -> bool {
                             if (value == 0 && rK0Short == region)
                               return true;
                             if (value == 1 && rLambda == region)
                               return true;
                             if (value == 2 && rAntiLambda == region)
                               return true;
                             return false;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassRegion, invMassRegion,
                           [](int rK0Short, int rLambda, int rAntiLambda, int value) -> int {
                             if (value == 0)
                               return rK0Short;
                             if (value == 1)
                               return rLambda;
                             if (value == 2)
                               return rAntiLambda;
                             return -1;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MCTrue, mcTrue,
                           [](int mcTrueK0Short, int mcTrueLambda, int mcTrueAntiLambda, int value) -> bool {
                             if (value == 0 && mcTrueK0Short == true)
                               return true;
                             if (value == 1 && mcTrueLambda == true)
                               return true;
                             if (value == 2 && mcTrueAntiLambda == true)
                               return true;
                             return false;
                           });
} // namespace assocV0s
DECLARE_SOA_TABLE(AssocV0s, "AOD", "ASSOCV0S", o2::soa::Index<>,
                  assocV0s::CollisionId, assocV0s::V0DataId,
                  assocV0s::CompatibleK0Short,
                  assocV0s::CompatibleLambda,
                  assocV0s::CompatibleAntiLambda,
                  assocV0s::MCTrueK0Short,
                  assocV0s::MCTrueLambda,
                  assocV0s::MCTrueAntiLambda,
                  assocV0s::MassRegionK0Short,
                  assocV0s::MassRegionLambda,
                  assocV0s::MassRegionAntiLambda,
                  assocV0s::Compatible<assocV0s::CompatibleK0Short, assocV0s::CompatibleLambda, assocV0s::CompatibleAntiLambda>,
                  assocV0s::InvMassRegionCheck<assocV0s::MassRegionK0Short, assocV0s::MassRegionLambda, assocV0s::MassRegionAntiLambda>,
                  assocV0s::InvMassRegion<assocV0s::MassRegionK0Short, assocV0s::MassRegionLambda, assocV0s::MassRegionAntiLambda>,
                  assocV0s::MCTrue<assocV0s::MCTrueK0Short, assocV0s::MCTrueLambda, assocV0s::MCTrueAntiLambda>);
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
DECLARE_SOA_COLUMN(MCTrueXiMinus, mcTrueXiMinus, bool);               // true XiMinus in mc
DECLARE_SOA_COLUMN(MCTrueXiPlus, mcTrueXiPlus, bool);                 // true XiPlus in mc
DECLARE_SOA_COLUMN(MCTrueOmegaMinus, mcTrueOmegaMinus, bool);         // true OmegaMinus in mc
DECLARE_SOA_COLUMN(MCTrueOmegaPlus, mcTrueOmegaPlus, bool);           // true OmegaPlus in mc
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
DECLARE_SOA_DYNAMIC_COLUMN(InvMassRegionCheck, invMassRegionCheck,
                           [](int rXi, int rOmega, int value, int region) -> bool {
                             if (value == 0 && rXi == region)
                               return true;
                             if (value == 1 && rXi == region)
                               return true;
                             if (value == 2 && rOmega == region)
                               return true;
                             if (value == 3 && rOmega == region)
                               return true;
                             return false;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(InvMassRegion, invMassRegion,
                           [](int rXi, int rOmega, int value) -> int {
                             if (value == 0 || value == 1)
                               return rXi;
                             if (value == 2 || value == 3)
                               return rOmega;
                             return -1;
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MCTrue, mcTrue,
                           [](int mcTrueXiMinus, int mcTrueXiPlus, int mcTrueOmegaMinus, int mcTrueOmegaPlus, int value) -> bool {
                             if (value == 0 && mcTrueXiMinus == true)
                               return true;
                             if (value == 1 && mcTrueXiPlus == true)
                               return true;
                             if (value == 2 && mcTrueOmegaMinus == true)
                               return true;
                             if (value == 3 && mcTrueOmegaPlus == true)
                               return true;
                             return false;
                           });
} // namespace assocCascades
DECLARE_SOA_TABLE(AssocCascades, "AOD", "ASSOCCASCADES", o2::soa::Index<>, assocCascades::CollisionId, assocCascades::CascDataId,
                  assocCascades::CompatibleXiMinus,
                  assocCascades::CompatibleXiPlus,
                  assocCascades::CompatibleOmegaMinus,
                  assocCascades::CompatibleOmegaPlus,
                  assocCascades::MCTrueXiMinus,
                  assocCascades::MCTrueXiPlus,
                  assocCascades::MCTrueOmegaMinus,
                  assocCascades::MCTrueOmegaPlus,
                  assocCascades::MassRegionXi,
                  assocCascades::MassRegionOmega,
                  assocCascades::Compatible<assocCascades::CompatibleXiMinus, assocCascades::CompatibleXiPlus, assocCascades::CompatibleOmegaMinus, assocCascades::CompatibleOmegaPlus>,
                  assocCascades::InvMassRegionCheck<assocCascades::MassRegionXi, assocCascades::MassRegionOmega>,
                  assocCascades::InvMassRegion<assocCascades::MassRegionXi, assocCascades::MassRegionOmega>,
                  assocCascades::MCTrue<assocCascades::MCTrueXiMinus, assocCascades::MCTrueXiPlus, assocCascades::MCTrueOmegaMinus, assocCascades::MCTrueOmegaPlus>);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHSTRANGECORRELATIONTABLES_H_
