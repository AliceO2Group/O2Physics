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
// delta-times
DECLARE_SOA_COLUMN(PosDeltaTimeLambdaPi, posDeltaTimeLambdaPi, float);   //! positive track DeltaTime from pion <- lambda expectation
DECLARE_SOA_COLUMN(PosDeltaTimeLambdaPr, posDeltaTimeLambdaPr, float);   //! positive track DeltaTime from proton <- lambda expectation
DECLARE_SOA_COLUMN(NegDeltaTimeLambdaPi, negDeltaTimeLambdaPi, float);   //! negative track DeltaTime from pion <- lambda expectation
DECLARE_SOA_COLUMN(NegDeltaTimeLambdaPr, negDeltaTimeLambdaPr, float);   //! negative track DeltaTime from proton <- lambda expectation
DECLARE_SOA_COLUMN(PosDeltaTimeK0ShortPi, posDeltaTimeK0ShortPi, float); //! positive track DeltaTime from pion <- k0short expectation
DECLARE_SOA_COLUMN(NegDeltaTimeK0ShortPi, negDeltaTimeK0ShortPi, float); //! positive track DeltaTime from pion <- k0short expectation

// n-sigmas
DECLARE_SOA_COLUMN(PosNSigmaLambdaPi, posNSigmaLambdaPi, float);   //! positive track NSigma from pion <- lambda expectation
DECLARE_SOA_COLUMN(PosNSigmaLambdaPr, posNSigmaLambdaPr, float);   //! positive track NSigma from proton <- lambda expectation
DECLARE_SOA_COLUMN(NegNSigmaLambdaPi, negNSigmaLambdaPi, float);   //! negative track NSigma from pion <- lambda expectation
DECLARE_SOA_COLUMN(NegNSigmaLambdaPr, negNSigmaLambdaPr, float);   //! negative track NSigma from proton <- lambda expectation
DECLARE_SOA_COLUMN(PosNSigmaK0ShortPi, posNSigmaK0ShortPi, float); //! positive track NSigma from pion <- k0short expectation
DECLARE_SOA_COLUMN(NegNSigmaK0ShortPi, negNSigmaK0ShortPi, float); //! positive track NSigma from pion <- k0short expectation
} // namespace v0data
DECLARE_SOA_TABLE(LaDeltaTimeTOF, "AOD", "LADELTATIMETOF",
                  v0data::PosDeltaTimeLambdaPi, v0data::PosDeltaTimeLambdaPr,
                  v0data::NegDeltaTimeLambdaPi, v0data::NegDeltaTimeLambdaPr);
DECLARE_SOA_TABLE(K0DeltaTimeTOF, "AOD", "K0DELTATIMETOF",
                  v0data::PosDeltaTimeK0ShortPi, v0data::NegDeltaTimeK0ShortPi);
DECLARE_SOA_TABLE(LaNSigmaTOF, "AOD", "LANSIGMATOF",
                  v0data::PosNSigmaLambdaPi, v0data::PosNSigmaLambdaPr,
                  v0data::NegNSigmaLambdaPi, v0data::NegNSigmaLambdaPr);
DECLARE_SOA_TABLE(K0NSigmaTOF, "AOD", "K0NSIGMATOF",
                  v0data::PosNSigmaK0ShortPi, v0data::NegNSigmaK0ShortPi);

namespace cascdata
{
// delta-times
DECLARE_SOA_COLUMN(PosDeltaTimeXiPi, posDeltaTimeXiPi, float);         //! positive track DeltaTime from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(PosDeltaTimeXiPr, posDeltaTimeXiPr, float);         //! positive track DeltaTime from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegDeltaTimeXiPi, negDeltaTimeXiPi, float);         //! negative track DeltaTime from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegDeltaTimeXiPr, negDeltaTimeXiPr, float);         //! negative track DeltaTime from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(BachDeltaTimeXiPi, bachDeltaTimeXiPi, float);       //! bachelor track DeltaTime from pion <- xi expectation
DECLARE_SOA_COLUMN(PosDeltaTimeOmegaPi, posDeltaTimeOmegaPi, float);   //! positive track DeltaTime from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(PosDeltaTimeOmegaPr, posDeltaTimeOmegaPr, float);   //! positive track DeltaTime from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegDeltaTimeOmegaPi, negDeltaTimeOmegaPi, float);   //! negative track DeltaTime from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegDeltaTimeOmegaPr, negDeltaTimeOmegaPr, float);   //! negative track DeltaTime from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(BachDeltaTimeOmegaPi, bachDeltaTimeOmegaPi, float); //! bachelor track DeltaTime from pion <- omega expectation

// n-sigmas
DECLARE_SOA_COLUMN(PosNSigmaXiPi, posNSigmaXiPi, float);         //! positive track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(PosNSigmaXiPr, posNSigmaXiPr, float);         //! positive track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegNSigmaXiPi, negNSigmaXiPi, float);         //! negative track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegNSigmaXiPr, negNSigmaXiPr, float);         //! negative track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(BachNSigmaXiPi, bachNSigmaXiPi, float);       //! bachelor track NSigma from pion <- xi expectation
DECLARE_SOA_COLUMN(PosNSigmaOmegaPi, posNSigmaOmegaPi, float);   //! positive track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(PosNSigmaOmegaPr, posNSigmaOmegaPr, float);   //! positive track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegNSigmaOmegaPi, negNSigmaOmegaPi, float);   //! negative track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegNSigmaOmegaPr, negNSigmaOmegaPr, float);   //! negative track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(BachNSigmaOmegaPi, bachNSigmaOmegaPi, float); //! bachelor track NSigma from pion <- omega expectation
} // namespace cascdata
DECLARE_SOA_TABLE(XiDeltaTimeTOF, "AOD", "XIDELTATIMETOF",
                  cascdata::PosDeltaTimeXiPi, cascdata::PosDeltaTimeXiPr,
                  cascdata::NegDeltaTimeXiPi, cascdata::NegDeltaTimeXiPr,
                  cascdata::BachDeltaTimeXiPi);
DECLARE_SOA_TABLE(OmDeltaTimeTOF, "AOD", "OMDELTATIMETOF",
                  cascdata::PosDeltaTimeOmegaPi, cascdata::PosDeltaTimeOmegaPr,
                  cascdata::NegDeltaTimeOmegaPi, cascdata::NegDeltaTimeOmegaPr,
                  cascdata::BachDeltaTimeOmegaPi);
DECLARE_SOA_TABLE(XiNSigmaTOF, "AOD", "XINSIGMATOF",
                  cascdata::PosNSigmaXiPi, cascdata::PosNSigmaXiPr,
                  cascdata::NegNSigmaXiPi, cascdata::NegNSigmaXiPr,
                  cascdata::BachNSigmaXiPi);
DECLARE_SOA_TABLE(OmNSigmaTOF, "AOD", "OmNSIGMATOF",
                  cascdata::PosNSigmaOmegaPi, cascdata::PosNSigmaOmegaPr,
                  cascdata::NegNSigmaOmegaPi, cascdata::NegNSigmaOmegaPr,
                  cascdata::BachNSigmaOmegaPi);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
