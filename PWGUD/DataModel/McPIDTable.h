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
#ifndef O2_ANALYSIS_UDDIFFMCSCAN_H
#define O2_ANALYSIS_UDDIFFMCSCAN_H

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace datascan
{
DECLARE_SOA_COLUMN(PID, pid, int);   //! MCTruth PID
DECLARE_SOA_COLUMN(Pt, pt, float);   //! MCTruth particle pt
DECLARE_SOA_COLUMN(sEL, sel, float); //! TPC nSigma el
DECLARE_SOA_COLUMN(sMU, smu, float); //! TPC nSigma mu
DECLARE_SOA_COLUMN(sPI, spi, float); //! TPC nSigma pi
DECLARE_SOA_COLUMN(sKA, ska, float); //! TPC nSigma ka
DECLARE_SOA_COLUMN(sPR, spr, float); //! TPC nSigma pr
} // namespace datascan

DECLARE_SOA_TABLE(UDnSigmas, "AOD", "UDNSIGMAS", //! MCTruth of particle PID and pt and TPC nSigma for el, mu, pi, ka, pr
                  datascan::PID, datascan::Pt,
                  datascan::sEL, datascan::sMU, datascan::sPI,
                  datascan::sKA, datascan::sPR);

} // namespace o2::aod

#endif // O2_ANALYSIS_UDDIFFMCSCAN_H
