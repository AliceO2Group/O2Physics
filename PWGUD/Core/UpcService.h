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

/// \file upcTrueGapService.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Apr/7/2026
/// \brief header defines some service outputs

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

namespace o2::aod
{
namespace upcservice
{
DECLARE_SOA_COLUMN(Truegapside, truegapside, int);
} // namespace upcservice
DECLARE_SOA_TABLE(Truegapside, "AOD", "TRUEGAPSIDE", upcservice::Truegapside);
} // namespace o2::aod
