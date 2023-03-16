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
#ifndef TUTORIALS_SRC_INTERMEDIATETABLES_H_
#define TUTORIALS_SRC_INTERMEDIATETABLES_H_
#include "Framework/AnalysisDataModel.h"

/// Definition of output table to store the flat tree for re-weighting model
/// training
namespace o2::aod
{
namespace vars
{
DECLARE_SOA_COLUMN(InPt, inpt, bool);
DECLARE_SOA_COLUMN(InEta, ineta, bool);
} // namespace vars

DECLARE_SOA_TABLE(Decisions, "AOD", "DECISIONS", o2::soa::Index<>,
                  vars::InPt,
                  vars::InEta);
} // namespace o2::aod

#endif // TUTORIALS_SRC_INTERMEDIATETABLES_H_
