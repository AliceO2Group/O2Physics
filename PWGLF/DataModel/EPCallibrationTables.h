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

/// \file EPCallibrationTables.h
///
/// \author Sourav Kundu <sourav.kundu@cern.ch>

#ifndef PWGLF_DATAMODEL_EPCALLIBRATIONTABLES_H_
#define PWGLF_DATAMODEL_EPCALLIBRATIONTABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace epcallibrationtable
{
DECLARE_SOA_COLUMN(TriggerEventEP, triggereventep, bool);
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(PsiFT0C, psiFT0C, float);
DECLARE_SOA_COLUMN(PsiFT0A, psiFT0A, float);
DECLARE_SOA_COLUMN(PsiTPC, psiTPC, float);
DECLARE_SOA_COLUMN(PsiTPCL, psiTPCL, float);
DECLARE_SOA_COLUMN(PsiTPCR, psiTPCR, float);
} // namespace epcallibrationtable
DECLARE_SOA_TABLE(EPCallibrationTables, "AOD", "EPCALLCOLS",
                  epcallibrationtable::TriggerEventEP,
                  epcallibrationtable::Cent,
                  epcallibrationtable::PsiFT0C,
                  epcallibrationtable::PsiFT0A,
                  epcallibrationtable::PsiTPC,
                  epcallibrationtable::PsiTPCL,
                  epcallibrationtable::PsiTPCR);
using EPCallibrationTable = EPCallibrationTables::iterator;
} // namespace o2::aod
#endif // PWGLF_DATAMODEL_EPCALLIBRATIONTABLES_H_
