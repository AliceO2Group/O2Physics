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
namespace dautrack
{
// ==== TPC INFORMATION ===
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float); //! track signal

DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float); //! Nsigma proton
DECLARE_SOA_COLUMN(TPCNSigmaHe, tpcNSigmaHe, float); //! Nsigma proton
} // namespace dautrack

DECLARE_SOA_TABLE(DauTrackTPCPIDs, "AOD", "DAUTRACKTPCPID", // nsigma table (for analysis)
                  dautrack::TPCSignal, dautrack::TPCNSigmaEl,
                  dautrack::TPCNSigmaPi, dautrack::TPCNSigmaKa,
                  dautrack::TPCNSigmaPr, dautrack::TPCNSigmaHe);

namespace v0data
{
// ==== TOF INFORMATION ===
// lengths as stored in the AO2D for TOF calculations
DECLARE_SOA_COLUMN(PosTOFLength, posTOFLength, float); //! positive track length
DECLARE_SOA_COLUMN(NegTOFLength, negTOFLength, float); //! negative track length

// delta-times
DECLARE_SOA_COLUMN(PosTOFDeltaTLaPi, posTOFDeltaTLaPi, float); //! positive track TOFDeltaT from pion <- lambda expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTLaPr, posTOFDeltaTLaPr, float); //! positive track TOFDeltaT from proton <- lambda expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTLaPi, negTOFDeltaTLaPi, float); //! negative track TOFDeltaT from pion <- lambda expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTLaPr, negTOFDeltaTLaPr, float); //! negative track TOFDeltaT from proton <- lambda expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTK0Pi, posTOFDeltaTK0Pi, float); //! positive track TOFDeltaT from pion <- k0short expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTK0Pi, negTOFDeltaTK0Pi, float); //! positive track TOFDeltaT from pion <- k0short expectation

// n-sigmas
DECLARE_SOA_COLUMN(PosNSigmaLaPi, posNSigmaLaPi, float); //! positive track NSigma from pion <- lambda expectation
DECLARE_SOA_COLUMN(PosNSigmaLaPr, posNSigmaLaPr, float); //! positive track NSigma from proton <- lambda expectation
DECLARE_SOA_COLUMN(NegNSigmaLaPi, negNSigmaLaPi, float); //! negative track NSigma from pion <- lambda expectation
DECLARE_SOA_COLUMN(NegNSigmaLaPr, negNSigmaLaPr, float); //! negative track NSigma from proton <- lambda expectation
DECLARE_SOA_COLUMN(PosNSigmaK0Pi, posNSigmaK0Pi, float); //! positive track NSigma from pion <- k0short expectation
DECLARE_SOA_COLUMN(NegNSigmaK0Pi, negNSigmaK0Pi, float); //! positive track NSigma from pion <- k0short expectation
} // namespace v0data

DECLARE_SOA_TABLE(V0TOF, "AOD", "V0TOF", // raw information table (for debug, etc)
                  v0data::PosTOFLength, v0data::NegTOFLength,
                  v0data::PosTOFDeltaTLaPi, v0data::PosTOFDeltaTLaPr,
                  v0data::NegTOFDeltaTLaPi, v0data::NegTOFDeltaTLaPr,
                  v0data::PosTOFDeltaTK0Pi, v0data::NegTOFDeltaTK0Pi);
DECLARE_SOA_TABLE(V0TOFPID, "AOD", "V0TOFPID", // nsigma table (for analysis)
                  v0data::PosNSigmaLaPi, v0data::PosNSigmaLaPr,
                  v0data::NegNSigmaLaPi, v0data::NegNSigmaLaPr,
                  v0data::PosNSigmaK0Pi, v0data::NegNSigmaK0Pi);

namespace cascdata
{
// ==== TOF INFORMATION ===
// lengths as stored in the AO2D for TOF calculations
DECLARE_SOA_COLUMN(PosTOFLength, posTOFLength, float);   //! positive track length
DECLARE_SOA_COLUMN(NegTOFLength, negTOFLength, float);   //! negative track length
DECLARE_SOA_COLUMN(BachTOFLength, bachTOFLength, float); //! bachelor track length

// delta-times
DECLARE_SOA_COLUMN(PosTOFDeltaTXiPi, posTOFDeltaTXiPi, float);   //! positive track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTXiPr, posTOFDeltaTXiPr, float);   //! positive track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTXiPi, negTOFDeltaTXiPi, float);   //! negative track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTXiPr, negTOFDeltaTXiPr, float);   //! negative track TOFDeltaT from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(BachTOFDeltaTXiPi, bachTOFDeltaTXiPi, float); //! bachelor track TOFDeltaT from pion <- xi expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTOmPi, posTOFDeltaTOmPi, float);   //! positive track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(PosTOFDeltaTOmPr, posTOFDeltaTOmPr, float);   //! positive track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTOmPi, negTOFDeltaTOmPi, float);   //! negative track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegTOFDeltaTOmPr, negTOFDeltaTOmPr, float);   //! negative track TOFDeltaT from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(BachTOFDeltaTOmPi, bachTOFDeltaTOmPi, float); //! bachelor track TOFDeltaT from pion <- omega expectation

// n-sigmas
DECLARE_SOA_COLUMN(PosNSigmaXiPi, posNSigmaXiPi, float);         //! positive track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(PosNSigmaXiPr, posNSigmaXiPr, float);         //! positive track NSigma from proton <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegNSigmaXiPi, negNSigmaXiPi, float);         //! negative track NSigma from pion <- lambda <- xi expectation
DECLARE_SOA_COLUMN(NegNSigmaXiPr, negNSigmaXiPr, float);         //! negative track NSigma from proton <- lambda <- xi expectation
DECLARE_SOA_COLUMN(BachNSigmaXiPi, bachNSigmaXiPi, float);       //! bachelor track NSigma from pion <- xi expectation
DECLARE_SOA_COLUMN(PosNSigmaOmPi, posNSigmaOmPi, float);         //! positive track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(PosNSigmaOmPr, posNSigmaOmPr, float);         //! positive track NSigma from proton <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegNSigmaOmPi, negNSigmaOmPi, float);         //! negative track NSigma from pion <- lambda <- omega expectation
DECLARE_SOA_COLUMN(NegNSigmaOmPr, negNSigmaOmPr, float);         //! negative track NSigma from proton <- lambda <- omega expectation
DECLARE_SOA_COLUMN(BachNSigmaOmKa, bachNSigmaOmKa, float);       //! bachelor track NSigma from kaon <- omega expectation
} // namespace cascdata

DECLARE_SOA_TABLE(CascTOF, "AOD", "CascTOF", // raw information table (for debug, etc)
                  cascdata::PosTOFLength, cascdata::NegTOFLength, cascdata::BachTOFLength,
                  cascdata::PosTOFDeltaTXiPi, cascdata::PosTOFDeltaTXiPr,
                  cascdata::NegTOFDeltaTXiPi, cascdata::NegTOFDeltaTXiPr,
                  cascdata::BachTOFDeltaTXiPi,
                  cascdata::PosTOFDeltaTOmPi, cascdata::PosTOFDeltaTOmPr,
                  cascdata::NegTOFDeltaTOmPi, cascdata::NegTOFDeltaTOmPr,
                  cascdata::BachTOFDeltaTOmPi);
DECLARE_SOA_TABLE(CascTOFPID, "AOD", "CASCTOFPID", // nsigma table (for analysis)
                  cascdata::PosNSigmaXiPi, cascdata::PosNSigmaXiPr,
                  cascdata::NegNSigmaXiPi, cascdata::NegNSigmaXiPr,
                  cascdata::BachNSigmaXiPi,
                  cascdata::PosNSigmaOmPi, cascdata::PosNSigmaOmPr,
                  cascdata::NegNSigmaOmPi, cascdata::NegNSigmaOmPr,
                  cascdata::BachNSigmaOmKa);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFSTRANGENESSPIDTABLES_H_
