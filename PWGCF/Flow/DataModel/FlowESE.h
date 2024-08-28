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

#ifndef PWGCF_FLOW_DATAMODEL_FLOWESE_H_
#define PWGCF_FLOW_DATAMODEL_FLOWESE_H_

#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

namespace o2::aod
{
namespace qVector
{
DECLARE_SOA_COLUMN(QFV0ARe, qFV0ARe, float);
DECLARE_SOA_COLUMN(QFV0AIm, qFV0AIm, float);
DECLARE_SOA_COLUMN(QFT0CRe, qFT0CRe, std::vector<float>);
DECLARE_SOA_COLUMN(QFT0CIm, qFT0CIm, std::vector<float>);
DECLARE_SOA_COLUMN(QPERCFT0C, qPERCFT0C, std::vector<float>);
DECLARE_SOA_COLUMN(QTEST, qTEST, float);
} // namespace qVector
DECLARE_SOA_TABLE(qVecFV0As, "AOD", "QVECFV0A", qVector::QFV0ARe, qVector::QFV0AIm);
DECLARE_SOA_TABLE(qVecFT0Cs, "AOD", "QVECFT0C", qVector::QFT0CRe, qVector::QFT0CIm);
DECLARE_SOA_TABLE(qPercentileFT0Cs, "AOD", "QPERCENTILEFT0C", qVector::QPERCFT0C);
using qVecFV0A = qVecFV0As::iterator;
using qVecFT0C = qVecFT0Cs::iterator;
using qPercentileFT0C = qPercentileFT0Cs::iterator;

} // namespace o2::aod

#endif // PWGCF_FLOW_DATAMODEL_FLOWESE_H_
