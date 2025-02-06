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
/// \brief Table definitions for selectors for writing out reduced data model for jets
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATASELECTOR_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATASELECTOR_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace jetreduceddataselector
{
DECLARE_SOA_COLUMN(IsCollisionSelected, isCollisionSelected, bool);
DECLARE_SOA_COLUMN(IsMcCollisionSelected, isMcCollisionSelected, bool);

} // namespace jetreduceddataselector
DECLARE_SOA_TABLE(JCollisionSelections, "AOD", "JCOLLSELECTION",
                  jetreduceddataselector::IsCollisionSelected);

DECLARE_SOA_TABLE(JMcCollisionSelections, "AOD", "JMCCOLLSELECTION",
                  jetreduceddataselector::IsMcCollisionSelected);
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATASELECTOR_H_
