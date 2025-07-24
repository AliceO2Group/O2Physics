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

#ifndef PWGMM_MULT_DATAMODEL_GENCENTRALITIES_H_
#define PWGMM_MULT_DATAMODEL_GENCENTRALITIES_H_

#include "Framework/AnalysisDataModel.h"
namespace o2::aod
{
namespace gencents
{
DECLARE_SOA_COLUMN(GenCentFT0C, gencentFT0C, float);
DECLARE_SOA_COLUMN(GenCentFT0M, gencentFT0M, float);
} // namespace gencents
DECLARE_SOA_TABLE(GenCents, "AOD", "GENCENT",
                  gencents::GenCentFT0C,
                  gencents::GenCentFT0M);
} // namespace o2::aod
#endif // PWGMM_MULT_DATAMODEL_GENCENTRALITIES_H_
