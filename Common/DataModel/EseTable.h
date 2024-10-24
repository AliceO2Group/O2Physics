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

// q vector framework with ESE (20/08/2024)
//
/// \author Joachim Hansen <joachim.hansen@cern.ch>
//

#ifndef COMMON_DATAMODEL_ESETABLE_H_
#define COMMON_DATAMODEL_ESETABLE_H_

#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

namespace o2::aod
{
namespace q_vector
{
DECLARE_SOA_COLUMN(QPERCFT0C, qPERCFT0C, std::vector<float>);
DECLARE_SOA_COLUMN(qPERCFT0A, qPERCFT0A, std::vector<float>);
DECLARE_SOA_COLUMN(qPERCFV0, qPERCFV0, std::vector<float>);
DECLARE_SOA_COLUMN(qPERCTPC, qPERCTPC, std::vector<float>);
DECLARE_SOA_COLUMN(FESECOL, fESECOL, std::vector<int>);
} // namespace q_vector
DECLARE_SOA_TABLE(QPercentileFT0Cs, "AOD", "QPERCENTILEFT0C", q_vector::QPERCFT0C);
DECLARE_SOA_TABLE(QPercentileFT0As, "AOD", "QPERCENTILEFT0A", q_vector::QPERCFT0A);
DECLARE_SOA_TABLE(QPercentileFV0s, "AOD", "QPERCENTILEFV0", q_vector::QPERCFV0);
DECLARE_SOA_TABLE(QPercentileTPCs, "AOD", "QPERCENTILETPC", q_vector::QPERCTPC);
DECLARE_SOA_TABLE(FEseCols, "AOD", "FEVENTSHAPE", q_vector::FESECOL);

using QPercentileFT0C = QPercentileFT0Cs::iterator;
using QPercentileFT0A = QPercentileFT0As::iterator;
using QPercentileFV0 = QPercentileFV0s::iterator;
using QPercentileTPC = QPercentileTPCs::iterator;
using FEseCol = FEseCols::iterator;

} // namespace o2::aod

#endif // COMMON_DATAMODEL_ESETABLE_H_
