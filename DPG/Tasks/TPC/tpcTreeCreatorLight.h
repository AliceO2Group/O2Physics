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

/// \file tpcTreeCreatorLight.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>

#ifndef DPG_TASKS_TPC_TPCTREECREATORLIGHT_H_
#define DPG_TASKS_TPC_TPCTREECREATORLIGHT_H_

#include "Common/Core/trackUtilities.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

enum ParticleSpecies {
  kPionTrack = BIT(0),
  kKaonTrack = BIT(1),
  kProtonTrack = BIT(2),
  kElectronTrack = BIT(3)
};

namespace o2::aod
{
namespace tpctree
{
DECLARE_SOA_COLUMN(NormMultTPC, normMultTPC, float);
DECLARE_SOA_COLUMN(NormNClustersTPC, normNClustersTPC, float);
DECLARE_SOA_COLUMN(PIDTrack, pidTrack, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace tpctree

DECLARE_SOA_TABLE(TPCTOFTree, "AOD", "TPCTOFTREE",
                  o2::aod::track::TPCSignal,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Tgl,
                  o2::aod::track::Signed1Pt,
                  o2::aod::track::Eta,
                  o2::aod::track::Phi,
                  o2::aod::track::Y,
                  tpctree::NormMultTPC,
                  tpctree::NormNClustersTPC,
                  tpctree::PIDTrack,
                  tpctree::RunNumber);
} // namespace o2::aod
#endif // DPG_TASKS_TPC_TPCTREECREATORLIGHT_H_
