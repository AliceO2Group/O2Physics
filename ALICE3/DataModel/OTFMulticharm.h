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

#ifndef ALICE3_DATAMODEL_OTFMULTICHARM_H_
#define ALICE3_DATAMODEL_OTFMULTICHARM_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otfmulticharm
{
DECLARE_SOA_INDEX_COLUMN_FULL(Cascade, cascade, int, UpgradeCascades, "_Cascade");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCPion1, xiCPion1, int, Tracks, "_Pi1XiC");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCPion2, xiCPion2, int, Tracks, "_Pi2XiC");
DECLARE_SOA_INDEX_COLUMN_FULL(XiCCPion, xiCCPion, int, Tracks, "_PiXiCC");

// topo vars
DECLARE_SOA_COLUMN(DCAXiCDaughters, dcaXiCDaughters, float);
DECLARE_SOA_COLUMN(DCAXiCCDaughters, dcaXiCCDaughters, float);

DECLARE_SOA_COLUMN(MXiC, mXiC, float);
DECLARE_SOA_COLUMN(MXiCC, mXiCC, float);

} // namespace otfmulticharm
DECLARE_SOA_TABLE(MultiCharmStates, "AOD", "MultiCharmStates",
                  o2::soa::Index<>,
                  otfcascade::CascadeId,
                  otfcascade::XiCPion1Id,
                  otfcascade::XiCPion2Id,
                  otfcascade::XiCCPionId,
                  otfcascade::DCAXiCDaughters,
                  otfcascade::DCAXiCCDaughters,
                  otfcascade::MXiC,
                  otfcascade::MXiCC);

using MultiCharmState = MultiCharmState::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFMULTICHARM_H_
