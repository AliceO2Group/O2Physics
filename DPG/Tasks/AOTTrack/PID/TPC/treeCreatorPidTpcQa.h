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

/// \file treeCreatorPidTpcQa.h
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Creates trees with PID QA variables along with variables used for NN training

#ifndef DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_
#define DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_

#include "DPG/Tasks/TPC/tpcSkimsTableCreator.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
DECLARE_SOA_TABLE(QaPidTpc, "AOD", "QAPIDTPC",
                  tpcskims::PidIndex,
                  tpcskims::Ft0Occ,
                  tpcskims::HadronicRate,
                  tpcskims::NormMultTPC,
                  tpcskims::NormNClustersTPC,
                  o2::aod::track::Phi,
                  o2::aod::track::Tgl,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Y)
} // namespace o2::aod

#endif // DPG_TASKS_AOTTRACK_PID_TPC_TREECREATORPIDTPCQA_H_
