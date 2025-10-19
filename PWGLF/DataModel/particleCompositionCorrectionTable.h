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
/// \file   ParticleCompositionCorrectionTable.h
/// \brief  Table for scaling MC particle abundances to match measurements
/// \author Mario Krüger <mario.kruger@cern.ch>
///

#ifndef PWGLF_DATAMODEL_PARTICLECOMPOSITIONCORRECTIONTABLE_H_
#define PWGLF_DATAMODEL_PARTICLECOMPOSITIONCORRECTIONTABLE_H_

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace PCC
{
DECLARE_SOA_COLUMN(PccWeight, pccWeight, float);
DECLARE_SOA_COLUMN(PccWeightSysUp, pccWeightSysUp, float);
DECLARE_SOA_COLUMN(PccWeightSysDown, pccWeightSysDown, float);
} // namespace PCC
DECLARE_SOA_TABLE(ParticleCompositionCorrection, "AOD", "PARTICLECOMPOSITIONCORRECTION", PCC::PccWeight, PCC::PccWeightSysUp, PCC::PccWeightSysDown);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_PARTICLECOMPOSITIONCORRECTIONTABLE_H_
