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

/// \file LFF1Tables.h
///
/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_FILTERF1PROTONTABLES_H_
#define PWGLF_DATAMODEL_FILTERF1PROTONTABLES_H_

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace filtering
{
DECLARE_SOA_COLUMN(TriggerEventF1Proton, triggereventf1proton, bool); //! F1 - proton femto trigger event
}
DECLARE_SOA_TABLE(F1ProtonFilters, "AOD", "F1ProtonFilters",
                  filtering::TriggerEventF1Proton);
using F1ProtonFilter = F1ProtonFilters::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_FILTERF1PROTONTABLES_H_
