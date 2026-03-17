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

#include <Framework/AnalysisDataModel.h>

#include <vector>

namespace o2::aod
{
namespace qvec
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(IsCalibrated, isCalibrated, bool);

DECLARE_SOA_COLUMN(QvecRe, qvecRe, std::vector<float>);
DECLARE_SOA_COLUMN(QvecIm, qvecIm, std::vector<float>);
DECLARE_SOA_COLUMN(QvecAmp, qvecAmp, std::vector<float>);

DECLARE_SOA_COLUMN(QvecFT0CReVec, qvecFT0CReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0CImVec, qvecFT0CImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0AReVec, qvecFT0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0AImVec, qvecFT0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0MReVec, qvecFT0MReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFT0MImVec, qvecFT0MImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFV0AReVec, qvecFV0AReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecFV0AImVec, qvecFV0AImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCposReVec, qvecTPCposReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCposImVec, qvecTPCposImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCnegReVec, qvecTPCnegReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCnegImVec, qvecTPCnegImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCallReVec, qvecTPCallReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecTPCallImVec, qvecTPCallImVec, std::vector<float>);

DECLARE_SOA_COLUMN(QvecFT0CRe, qvecFT0CRe, float);
DECLARE_SOA_COLUMN(QvecFT0CIm, qvecFT0CIm, float);
DECLARE_SOA_COLUMN(QvecFT0ARe, qvecFT0ARe, float);
DECLARE_SOA_COLUMN(QvecFT0AIm, qvecFT0AIm, float);
DECLARE_SOA_COLUMN(QvecFT0MRe, qvecFT0MRe, float);
DECLARE_SOA_COLUMN(QvecFT0MIm, qvecFT0MIm, float);
DECLARE_SOA_COLUMN(QvecFV0ARe, qvecFV0ARe, float);
DECLARE_SOA_COLUMN(QvecFV0AIm, qvecFV0AIm, float);
DECLARE_SOA_COLUMN(QvecTPCposRe, qvecTPCposRe, float);
DECLARE_SOA_COLUMN(QvecTPCposIm, qvecTPCposIm, float);
DECLARE_SOA_COLUMN(QvecTPCnegRe, qvecTPCnegRe, float);
DECLARE_SOA_COLUMN(QvecTPCnegIm, qvecTPCnegIm, float);
DECLARE_SOA_COLUMN(QvecTPCallRe, qvecTPCallRe, float);
DECLARE_SOA_COLUMN(QvecTPCallIm, qvecTPCallIm, float);

DECLARE_SOA_COLUMN(SumAmplFT0C, sumAmplFT0C, float);
DECLARE_SOA_COLUMN(SumAmplFT0A, sumAmplFT0A, float);
DECLARE_SOA_COLUMN(SumAmplFT0M, sumAmplFT0M, float);
DECLARE_SOA_COLUMN(SumAmplFV0A, sumAmplFV0A, float);
DECLARE_SOA_COLUMN(NTrkTPCpos, nTrkTPCpos, int);
DECLARE_SOA_COLUMN(NTrkTPCneg, nTrkTPCneg, int);
DECLARE_SOA_COLUMN(NTrkTPCall, nTrkTPCall, int);
DECLARE_SOA_COLUMN(LabelsTPCpos, labelsTPCpos, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsTPCneg, labelsTPCneg, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsTPCall, labelsTPCall, std::vector<int>);

// Deprecated, will be removed in future after transition time //
DECLARE_SOA_COLUMN(QvecBPosReVec, qvecBPosReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBPosImVec, qvecBPosImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBNegReVec, qvecBNegReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBNegImVec, qvecBNegImVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBTotReVec, qvecBTotReVec, std::vector<float>);
DECLARE_SOA_COLUMN(QvecBTotImVec, qvecBTotImVec, std::vector<float>);

DECLARE_SOA_COLUMN(QvecBPosRe, qvecBPosRe, float);
DECLARE_SOA_COLUMN(QvecBPosIm, qvecBPosIm, float);
DECLARE_SOA_COLUMN(QvecBNegRe, qvecBNegRe, float);
DECLARE_SOA_COLUMN(QvecBNegIm, qvecBNegIm, float);
DECLARE_SOA_COLUMN(QvecBTotRe, qvecBTotRe, float);
DECLARE_SOA_COLUMN(QvecBTotIm, qvecBTotIm, float);

DECLARE_SOA_COLUMN(NTrkBPos, nTrkBPos, int);
DECLARE_SOA_COLUMN(NTrkBNeg, nTrkBNeg, int);
DECLARE_SOA_COLUMN(NTrkBTot, nTrkBTot, int);
DECLARE_SOA_COLUMN(LabelsBPos, labelsBPos, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsBNeg, labelsBNeg, std::vector<int>);
DECLARE_SOA_COLUMN(LabelsBTot, labelsBTot, std::vector<int>);
/////////////////////////////////////////////////////////////////
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORDEVS", //! Table with all Qvectors.
                  qvec::Cent, qvec::IsCalibrated, qvec::QvecRe, qvec::QvecIm, qvec::QvecAmp);
using Qvector = Qvectors::iterator;

DECLARE_SOA_TABLE(QvectorFT0Cs, "AOD", "QVECTORSFT0C", qvec::IsCalibrated, qvec::QvecFT0CRe, qvec::QvecFT0CIm, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(QvectorFT0As, "AOD", "QVECTORSFT0A", qvec::IsCalibrated, qvec::QvecFT0ARe, qvec::QvecFT0AIm, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(QvectorFT0Ms, "AOD", "QVECTORSFT0M", qvec::IsCalibrated, qvec::QvecFT0MRe, qvec::QvecFT0MIm, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(QvectorFV0As, "AOD", "QVECTORSFV0A", qvec::IsCalibrated, qvec::QvecFV0ARe, qvec::QvecFV0AIm, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(QvectorTPCposs, "AOD", "QVECTORSTPCPOS", qvec::IsCalibrated, qvec::QvecTPCposRe, qvec::QvecTPCposIm, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(QvectorTPCnegs, "AOD", "QVECTORSTPCNEG", qvec::IsCalibrated, qvec::QvecTPCnegRe, qvec::QvecTPCnegIm, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(QvectorTPCalls, "AOD", "QVECTORSTPCALL", qvec::IsCalibrated, qvec::QvecTPCallRe, qvec::QvecTPCallIm, qvec::NTrkTPCall, qvec::LabelsTPCall);

DECLARE_SOA_TABLE(QvectorFT0CVecs, "AOD", "QVECTORSFT0CVEC", qvec::IsCalibrated, qvec::QvecFT0CReVec, qvec::QvecFT0CImVec, qvec::SumAmplFT0C);
DECLARE_SOA_TABLE(QvectorFT0AVecs, "AOD", "QVECTORSFT0AVEC", qvec::IsCalibrated, qvec::QvecFT0AReVec, qvec::QvecFT0AImVec, qvec::SumAmplFT0A);
DECLARE_SOA_TABLE(QvectorFT0MVecs, "AOD", "QVECTORSFT0MVEC", qvec::IsCalibrated, qvec::QvecFT0MReVec, qvec::QvecFT0MImVec, qvec::SumAmplFT0M);
DECLARE_SOA_TABLE(QvectorFV0AVecs, "AOD", "QVECTORSFV0AVEC", qvec::IsCalibrated, qvec::QvecFV0AReVec, qvec::QvecFV0AImVec, qvec::SumAmplFV0A);
DECLARE_SOA_TABLE(QvectorTPCposVecs, "AOD", "QVECTORSTPCPVEC", qvec::IsCalibrated, qvec::QvecTPCposReVec, qvec::QvecTPCposImVec, qvec::NTrkTPCpos, qvec::LabelsTPCpos);
DECLARE_SOA_TABLE(QvectorTPCnegVecs, "AOD", "QVECTORSTPCNVEC", qvec::IsCalibrated, qvec::QvecTPCnegReVec, qvec::QvecTPCnegImVec, qvec::NTrkTPCneg, qvec::LabelsTPCneg);
DECLARE_SOA_TABLE(QvectorTPCallVecs, "AOD", "QVECTORSTPCAVEC", qvec::IsCalibrated, qvec::QvecTPCallReVec, qvec::QvecTPCallImVec, qvec::NTrkTPCall, qvec::LabelsTPCall);

using QvectorFT0C = QvectorFT0Cs::iterator;
using QvectorFT0A = QvectorFT0As::iterator;
using QvectorFT0M = QvectorFT0Ms::iterator;
using QvectorFV0A = QvectorFV0As::iterator;
using QvectorTPCpos = QvectorTPCposs::iterator;
using QvectorTPCneg = QvectorTPCnegs::iterator;
using QvectorTPCall = QvectorTPCalls::iterator;

using QvectorFT0CVec = QvectorFT0CVecs::iterator;
using QvectorFT0AVec = QvectorFT0AVecs::iterator;
using QvectorFT0MVec = QvectorFT0MVecs::iterator;
using QvectorFV0AVec = QvectorFV0AVecs::iterator;
using QvectorTPCposVec = QvectorTPCposVecs::iterator;
using QvectorTPCnegVec = QvectorTPCnegVecs::iterator;
using QvectorTPCallVec = QvectorTPCallVecs::iterator;

// Deprecated, will be removed in future after transition time //
DECLARE_SOA_TABLE(QvectorBPoss, "AOD", "QVECTORSBPOS", qvec::IsCalibrated, qvec::QvecBPosRe, qvec::QvecBPosIm, qvec::NTrkBPos, qvec::LabelsBPos);
DECLARE_SOA_TABLE(QvectorBNegs, "AOD", "QVECTORSBNEG", qvec::IsCalibrated, qvec::QvecBNegRe, qvec::QvecBNegIm, qvec::NTrkBNeg, qvec::LabelsBNeg);
DECLARE_SOA_TABLE(QvectorBTots, "AOD", "QVECTORSBTOT", qvec::IsCalibrated, qvec::QvecBTotRe, qvec::QvecBTotIm, qvec::NTrkBTot, qvec::LabelsBTot);

DECLARE_SOA_TABLE(QvectorBPosVecs, "AOD", "QVECTORSBPOSVEC", qvec::IsCalibrated, qvec::QvecBPosReVec, qvec::QvecBPosImVec, qvec::NTrkBPos, qvec::LabelsBPos);
DECLARE_SOA_TABLE(QvectorBNegVecs, "AOD", "QVECTORSBNEGVEC", qvec::IsCalibrated, qvec::QvecBNegReVec, qvec::QvecBNegImVec, qvec::NTrkBNeg, qvec::LabelsBNeg);
DECLARE_SOA_TABLE(QvectorBTotVecs, "AOD", "QVECTORSBTOTVEC", qvec::IsCalibrated, qvec::QvecBTotReVec, qvec::QvecBTotImVec, qvec::NTrkBTot, qvec::LabelsBTot);

using QvectorBPos = QvectorBPoss::iterator;
using QvectorBNeg = QvectorBNegs::iterator;
using QvectorBTot = QvectorBTots::iterator;

using QvectorBPosVec = QvectorBPosVecs::iterator;
using QvectorBNegVec = QvectorBNegVecs::iterator;
using QvectorBTotVec = QvectorBTotVecs::iterator;
/////////////////////////////////////////////////////////////////

} // namespace o2::aod

#endif // COMMON_DATAMODEL_QVECTORS_H_
