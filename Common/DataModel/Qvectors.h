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
/// \file   Qvectors.h
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Declaration of the table for the (un)corrected Q-vectors for the event plane
/// determination.
///

#ifndef COMMON_DATAMODEL_QVECTORS_H_
#define COMMON_DATAMODEL_QVECTORS_H_

#include <vector>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace qvec
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(CentBin, centBin, int);
DECLARE_SOA_COLUMN(QvecFT0CUncorRe, qvecFT0CUncorRe, float);
DECLARE_SOA_COLUMN(QvecFT0CUncorIm, qvecFT0CUncorIm, float);
DECLARE_SOA_COLUMN(QvecFT0CRectrRe, qvecFT0CRectrRe, float);
DECLARE_SOA_COLUMN(QvecFT0CRectrIm, qvecFT0CRectrIm, float);
DECLARE_SOA_COLUMN(QvecFT0CTwistRe, qvecFT0CTwistRe, float);
DECLARE_SOA_COLUMN(QvecFT0CTwistIm, qvecFT0CTwistIm, float);
DECLARE_SOA_COLUMN(QvecFT0CRe, qvecFT0CRe, float);
DECLARE_SOA_COLUMN(QvecFT0CIm, qvecFT0CIm, float);

DECLARE_SOA_COLUMN(QvecFT0AUncorRe, qvecFT0AUncorRe, float);
DECLARE_SOA_COLUMN(QvecFT0AUncorIm, qvecFT0AUncorIm, float);
DECLARE_SOA_COLUMN(QvecFT0ARectrRe, qvecFT0ARectrRe, float);
DECLARE_SOA_COLUMN(QvecFT0ARectrIm, qvecFT0ARectrIm, float);
DECLARE_SOA_COLUMN(QvecFT0ATwistRe, qvecFT0ATwistRe, float);
DECLARE_SOA_COLUMN(QvecFT0ATwistIm, qvecFT0ATwistIm, float);
DECLARE_SOA_COLUMN(QvecFT0ARe, qvecFT0ARe, float);
DECLARE_SOA_COLUMN(QvecFT0AIm, qvecFT0AIm, float);

DECLARE_SOA_COLUMN(QvecFT0MUncorRe, qvecFT0MUncorRe, float);
DECLARE_SOA_COLUMN(QvecFT0MUncorIm, qvecFT0MUncorIm, float);
DECLARE_SOA_COLUMN(QvecFT0MRectrRe, qvecFT0MRectrRe, float);
DECLARE_SOA_COLUMN(QvecFT0MRectrIm, qvecFT0MRectrIm, float);
DECLARE_SOA_COLUMN(QvecFT0MTwistRe, qvecFT0MTwistRe, float);
DECLARE_SOA_COLUMN(QvecFT0MTwistIm, qvecFT0MTwistIm, float);
DECLARE_SOA_COLUMN(QvecFT0MRe, qvecFT0MRe, float);
DECLARE_SOA_COLUMN(QvecFT0MIm, qvecFT0MIm, float);

DECLARE_SOA_COLUMN(QvecFV0AUncorRe, qvecFV0AUncorRe, float);
DECLARE_SOA_COLUMN(QvecFV0AUncorIm, qvecFV0AUncorIm, float);
DECLARE_SOA_COLUMN(QvecFV0ARectrRe, qvecFV0ARectrRe, float);
DECLARE_SOA_COLUMN(QvecFV0ARectrIm, qvecFV0ARectrIm, float);
DECLARE_SOA_COLUMN(QvecFV0ATwistRe, qvecFV0ATwistRe, float);
DECLARE_SOA_COLUMN(QvecFV0ATwistIm, qvecFV0ATwistIm, float);
DECLARE_SOA_COLUMN(QvecFV0ARe, qvecFV0ARe, float);
DECLARE_SOA_COLUMN(QvecFV0AIm, qvecFV0AIm, float);

DECLARE_SOA_COLUMN(QvecBPosUncorRe, qvecBPosUncorRe, float);
DECLARE_SOA_COLUMN(QvecBPosUncorIm, qvecBPosUncorIm, float);
DECLARE_SOA_COLUMN(QvecBPosRectrRe, qvecBPosRectrRe, float);
DECLARE_SOA_COLUMN(QvecBPosRectrIm, qvecBPosRectrIm, float);
DECLARE_SOA_COLUMN(QvecBPosTwistRe, qvecBPosTwistRe, float);
DECLARE_SOA_COLUMN(QvecBPosTwistIm, qvecBPosTwistIm, float);
DECLARE_SOA_COLUMN(QvecBPosRe, qvecBPosRe, float);
DECLARE_SOA_COLUMN(QvecBPosIm, qvecBPosIm, float);

DECLARE_SOA_COLUMN(QvecBNegUncorRe, qvecBNegUncorRe, float);
DECLARE_SOA_COLUMN(QvecBNegUncorIm, qvecBNegUncorIm, float);
DECLARE_SOA_COLUMN(QvecBNegRectrRe, qvecBNegRectrRe, float);
DECLARE_SOA_COLUMN(QvecBNegRectrIm, qvecBNegRectrIm, float);
DECLARE_SOA_COLUMN(QvecBNegTwistRe, qvecBNegTwistRe, float);
DECLARE_SOA_COLUMN(QvecBNegTwistIm, qvecBNegTwistIm, float);
DECLARE_SOA_COLUMN(QvecBNegRe, qvecBNegRe, float);
DECLARE_SOA_COLUMN(QvecBNegIm, qvecBNegIm, float);

DECLARE_SOA_COLUMN(SumAmplFT0C, sumAmplFT0C, float);
DECLARE_SOA_COLUMN(SumAmplFT0A, sumAmplFT0A, float);
DECLARE_SOA_COLUMN(SumAmplFT0M, sumAmplFT0M, float);
DECLARE_SOA_COLUMN(SumAmplFV0A, sumAmplFV0A, float);
DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
DECLARE_SOA_COLUMN(LabelsBPos, labelsBPos, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsBNeg, labelsBNeg, std::vector<int>);
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORDEVS", //! Table with all Qvectors.
                  qvec::Cent, qvec::CentBin,
                  qvec::QvecFT0CUncorRe, qvec::QvecFT0CUncorIm, qvec::QvecFT0CRectrRe, qvec::QvecFT0CRectrIm,
                  qvec::QvecFT0CTwistRe, qvec::QvecFT0CTwistIm, qvec::QvecFT0CRe, qvec::QvecFT0CIm,
                  qvec::QvecFT0AUncorRe, qvec::QvecFT0AUncorIm, qvec::QvecFT0ARectrRe, qvec::QvecFT0ARectrIm,
                  qvec::QvecFT0ATwistRe, qvec::QvecFT0ATwistIm, qvec::QvecFT0ARe, qvec::QvecFT0AIm,
                  qvec::QvecFT0MUncorRe, qvec::QvecFT0MUncorIm, qvec::QvecFT0MRectrRe, qvec::QvecFT0MRectrIm,
                  qvec::QvecFT0MTwistRe, qvec::QvecFT0MTwistIm, qvec::QvecFT0MRe, qvec::QvecFT0MIm,
                  qvec::QvecFV0AUncorRe, qvec::QvecFV0AUncorIm, qvec::QvecFV0ARectrRe, qvec::QvecFV0ARectrIm,
                  qvec::QvecFV0ATwistRe, qvec::QvecFV0ATwistIm, qvec::QvecFV0ARe, qvec::QvecFV0AIm,
                  qvec::QvecBPosUncorRe, qvec::QvecBPosUncorIm, qvec::QvecBPosRectrRe, qvec::QvecBPosRectrIm,
                  qvec::QvecBPosTwistRe, qvec::QvecBPosTwistIm, qvec::QvecBPosRe, qvec::QvecBPosIm,
                  qvec::QvecBNegUncorRe, qvec::QvecBNegUncorIm, qvec::QvecBNegRectrRe, qvec::QvecBNegRectrIm,
                  qvec::QvecBNegTwistRe, qvec::QvecBNegTwistIm, qvec::QvecBNegRe, qvec::QvecBNegIm,
                  qvec::SumAmplFT0C, qvec::SumAmplFT0A, qvec::SumAmplFT0M, qvec::SumAmplFV0A, qvec::NTrkBPos, qvec::NTrkBNeg,
                  qvec::LabelsBPos, qvec::LabelsBNeg);
using Qvector = Qvectors::iterator;

DECLARE_SOA_TABLE(QvectorFT0Cs, "AOD", "QVECTORSFT0C", qvec::CentBin, qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(QvectorFT0As, "AOD", "QVECTORSFT0A", qvec::CentBin, qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(QvectorFT0Ms, "AOD", "QVECTORSFT0M", qvec::CentBin, qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(QvectorFV0As, "AOD", "QVECTORSFV0A", qvec::CentBin, qvec::QvecFV0ARe, qvec::QvecFV0AIm, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(QvectorBPoss, "AOD", "QVECTORSBPOS", qvec::CentBin, qvec::QvecBPosRe, qvec::QvecBPosIm, qvec::NTrkBPos, qvec::LabelsBPos);
DECLARE_SOA_TABLE(QvectorBNegs, "AOD", "QVECTORSBNEG", qvec::CentBin, qvec::QvecBNegRe, qvec::QvecBNegIm, qvec::NTrkBNeg, qvec::LabelsBNeg);

using QvectorFT0C = QvectorFT0Cs::iterator;
using QvectorFT0A = QvectorFT0As::iterator;
using QvectorFT0M = QvectorFT0Ms::iterator;
using QvectorFV0A = QvectorFV0As::iterator;
using QvectorBPos = QvectorBPoss::iterator;
using QvectorBNeg = QvectorBNegs::iterator;

} // namespace o2::aod

#endif // COMMON_DATAMODEL_QVECTORS_H_
