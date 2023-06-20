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

// Defines TOF PID tables for strangeness.
// Entries calculated per candidate, tables are joinable with v0/cascdata tables.

#ifndef PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
#define PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::aod
{
namespace v0data
{
DECLARE_SOA_COLUMN(PosDeltaTimePi, posDeltaTimePi, float); //! positive track DeltaTime from pion expectation
DECLARE_SOA_COLUMN(PosDeltaTimePr, posDeltaTimePr, float); //! positive track DeltaTime from proton expectation
DECLARE_SOA_COLUMN(NegDeltaTimePi, negDeltaTimePi, float); //! negative track DeltaTime from pion expectation
DECLARE_SOA_COLUMN(NegDeltaTimePr, negDeltaTimePr, float); //! negative track DeltaTime from proton expectation
DECLARE_SOA_COLUMN(PosNSigmaPi, posNSigmaPi, float);       //! positive track Nsigma from pion expectation
DECLARE_SOA_COLUMN(PosNSigmaPr, posNSigmaPr, float);       //! positive track Nsigma from proton expectation
DECLARE_SOA_COLUMN(NegNSigmaPi, negNSigmaPi, float);       //! negative track Nsigma from pion expectation
DECLARE_SOA_COLUMN(NegNSigmaPr, negNSigmaPr, float);       //! negative track Nsigma from proton expectation
} // namespace v0data
DECLARE_SOA_TABLE(V0DeltaTimeTOF, "AOD", "V0DELTATIMETOF",
                  v0data::PosDeltaTimePi, v0data::PosDeltaTimePr,
                  v0data::NegDeltaTimePi, v0data::NegDeltaTimePr);
DECLARE_SOA_TABLE(V0NSigmaTOF, "AOD", "V0NSIGMATOF",
                  v0data::PosNSigmaPi, v0data::PosNSigmaPr,
                  v0data::NegNSigmaPi, v0data::NegNSigmaPr);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
